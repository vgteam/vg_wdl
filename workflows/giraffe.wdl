version 1.0

### giraffe.wdl ###
## Author: Charles Markello, Jean Monlong, Adam Novak
## Description: Core VG Giraffe mapping, usable for DeepVariant.
## Reference: https://github.com/vgteam/vg/wiki

import "../tasks/bioinfo_utils.wdl" as utils

workflow Giraffe {
    input {
        File? INPUT_READ_FILE_1                         # Input sample 1st read pair fastq.gz
        File? INPUT_READ_FILE_2                         # Input sample 2nd read pair fastq.gz
        File? INPUT_CRAM_FILE                           # Input CRAM file
        File? CRAM_REF                                  # Genome fasta file associated with the CRAM file
        File? CRAM_REF_INDEX                            # Index of the fasta file associated with the CRAM file
        String SAMPLE_NAME                              # The sample name
        Int MAX_FRAGMENT_LENGTH = 3000                  # Maximum distance at which to mark paired reads properly paired
        String VG_CONTAINER = "quay.io/vgteam/vg:v1.43.0" # VG Container used in the pipeline
        Int READS_PER_CHUNK = 20000000                  # Number of reads contained in each mapping chunk (20000000 for wgs)
        String GIRAFFE_OPTIONS = ""                     # (OPTIONAL) extra command line options for Giraffe mapper
        Array[String]+? CONTIGS                         # (OPTIONAL) Desired reference genome contigs, which are all paths in the XG index.
        File? PATH_LIST_FILE                            # (OPTIONAL) Text file where each line is a path name in the XG index, to use instead of CONTIGS. If neither is given, paths are extracted from the XG and subset to chromosome-looking paths.
        String REFERENCE_PREFIX = ""                    # Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)
        File XG_FILE                                    # Path to .xg index file
        File GBWT_FILE                                  # Path to .gbwt index file
        File GGBWT_FILE                                 # Path to .gg index file
        File DIST_FILE                                  # Path to .dist index file
        File MIN_FILE                                   # Path to .min index file
        Boolean LEFTALIGN_BAM = true                    # Whether or not to left-align reads in the BAM before DV
        Boolean REALIGN_INDELS = true                   # Whether or not to realign reads near indels before DV
        Int REALIGNMENT_EXPANSION_BASES = 160           # Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions.
        Int MIN_MAPQ = 1                                # Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong.
        Boolean OUTPUT_SINGLE_BAM = true                # Should a single merged BAM file be saved?
        Boolean OUTPUT_CALLING_BAMS = false             # Should individual contig BAMs be saved?
        Boolean OUTPUT_GAF = false                      # Should a GAF file with the aligned reads be saved?
        Boolean OUTPUT_GAM = true                       # Should a GAM file with the aligned reads be saved?
        Int SPLIT_READ_CORES = 8
        Int SPLIT_READ_DISK = 60
        Int MAP_CORES = 16
        Int MAP_DISK = 200
        Int MAP_MEM = 120
        Int POSTPROCESS_DISK = 40
        File? REFERENCE_FILE                            # (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference.
        File? REFERENCE_INDEX_FILE                      # (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.
        File? REFERENCE_DICT_FILE                       # (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if REFERENCE_INDEX_FILE is set. 
    }

    if(defined(INPUT_CRAM_FILE) && defined(CRAM_REF) && defined(CRAM_REF_INDEX)) {
	call convertCRAMtoFASTQ {
            input:
            in_cram_file=INPUT_CRAM_FILE,
            in_ref_file=CRAM_REF,
            in_ref_index_file=CRAM_REF_INDEX,
            in_cores=SPLIT_READ_CORES
	}
    }

    File read_1_file = select_first([INPUT_READ_FILE_1, convertCRAMtoFASTQ.output_fastq_1_file])
    File read_2_file = select_first([INPUT_READ_FILE_2, convertCRAMtoFASTQ.output_fastq_2_file])
    
    # Split input reads into chunks for parallelized mapping
    call splitReads as firstReadPair {
        input:
            in_read_file=read_1_file,
            in_pair_id="1",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES,
            in_split_read_disk=SPLIT_READ_DISK
    }
    call splitReads as secondReadPair {
        input:
            in_read_file=read_2_file,
            in_pair_id="2",
            in_reads_per_chunk=READS_PER_CHUNK,
            in_split_read_cores=SPLIT_READ_CORES,
            in_split_read_disk=SPLIT_READ_DISK
    }

    # Which path names to work on?
    if (!defined(CONTIGS)) {
        if (!defined(PATH_LIST_FILE)) {
            # Extract path names to call against from xg file if PATH_LIST_FILE input not provided
            # Filter down to major paths, because GRCh38 includes thousands of
            # decoys and unplaced/unlocalized contigs, and we can't efficiently
            # scatter across them, nor do we care about accuracy on them, and also
            # calling on the decoys is semantically meaningless.
            call extractSubsetPathNames {
                input:
                    in_xg_file=XG_FILE,
                    in_vg_container=VG_CONTAINER,
                    in_extract_disk=MAP_DISK,
                    in_extract_mem=MAP_MEM
            }
        }
    } 
    if (defined(CONTIGS)) {
        # Put the paths in a file to use later. We know the value is defined,
        # but WDL is a bit low on unboxing calls for optionals so we use
        # select_first.
        File written_path_names_file = write_lines(select_first([CONTIGS]))
    }
    File pipeline_path_list_file = select_first([PATH_LIST_FILE, extractSubsetPathNames.output_path_list_file, written_path_names_file])
    
    # To make sure that we have a FASTA reference with a contig set that
    # exactly matches the graph, we generate it ourselves, from the graph.
    if (!defined(REFERENCE_FILE)) {
        call extractReference {
            input:
                in_xg_file=XG_FILE,
                in_path_list_file=pipeline_path_list_file,
                in_vg_container=VG_CONTAINER,
                in_extract_disk=MAP_DISK,
                in_extract_mem=MAP_MEM
        }
    }
    File reference_file = select_first([REFERENCE_FILE, extractReference.reference_file])
    
    if (!defined(REFERENCE_INDEX_FILE)) {
        call utils.indexReference {
            input:
                in_reference_file=reference_file,
                in_index_disk=MAP_DISK,
                in_index_mem=MAP_MEM
        }
    }
    File reference_index_file = select_first([REFERENCE_INDEX_FILE, indexReference.reference_index_file])
    File reference_dict_file = select_first([REFERENCE_DICT_FILE, indexReference.reference_dict_file])

    ################################################################
    # Distribute vg mapping operation over each chunked read pair #
    ################################################################
    Array[Pair[File,File]] read_pair_chunk_files_list = zip(firstReadPair.output_read_chunks, secondReadPair.output_read_chunks)
    scatter (read_pair_chunk_files in read_pair_chunk_files_list) {
        call runVGGIRAFFE {
            input:
                in_left_read_pair_chunk_file=read_pair_chunk_files.left,
                in_right_read_pair_chunk_file=read_pair_chunk_files.right,
                in_vg_container=VG_CONTAINER,
                in_giraffe_options=GIRAFFE_OPTIONS,
                in_xg_file=XG_FILE,
                in_gbwt_file=GBWT_FILE,
                in_ggbwt_file=GGBWT_FILE,
                in_dist_file=DIST_FILE,
                in_min_file=MIN_FILE,
                # We always need to pass a full dict file here, with lengths,
                # because if we pass just path lists and the paths are not
                # completely contained in the graph (like if we're working on
                # GRCh38 paths in a CHM13-based graph), giraffe won't be able
                # to get the path lengths and will crash.
                # TODO: Somehow this problem is supposed to go away if we pull
                # any GRCh38. prefix off the path names by setting
                # REFERENCE_PREFIX and making sure the prefix isn't in the
                # truth set.
                # See <https://github.com/adamnovak/giraffe-dv-wdl/pull/2#issuecomment-955096920>
                in_sample_name=SAMPLE_NAME,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
        call surjectGAMtoBAM {
                input:
                in_gam_file=runVGGIRAFFE.chunk_gam_file,
                in_xg_file=XG_FILE,
                in_path_list_file=pipeline_path_list_file,
                in_sample_name=SAMPLE_NAME,
                in_max_fragment_length=MAX_FRAGMENT_LENGTH,
                in_vg_container=VG_CONTAINER,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
            }
        if (REFERENCE_PREFIX != "") {
            # use samtools to replace the header contigs with those from our dict.
            # this allows the header to contain contigs that are not in the graph,
            # which is more general and lets CHM13-based graphs be used to call on GRCh38
            # also, strip out contig prefixes in the BAM body
            call fixBAMContigNaming {
                input:
                    in_bam_file=surjectGAMtoBAM.chunk_bam_file,
                    in_ref_dict=reference_dict_file,
                    in_prefix_to_strip=REFERENCE_PREFIX,
                    in_map_cores=MAP_CORES,
                    in_map_disk=MAP_DISK,
                    in_map_mem=MAP_MEM
            }
        }
        File properly_named_bam_file = select_first([fixBAMContigNaming.fixed_bam_file, surjectGAMtoBAM.chunk_bam_file]) 
        call sortBAMFile {
            input:
                in_sample_name=SAMPLE_NAME,
                in_bam_chunk_file=properly_named_bam_file,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
    }
    Array[File] alignment_chunk_bam_files = select_all(sortBAMFile.sorted_chunk_bam)

    call mergeAlignmentBAMChunks {
        input:
            in_sample_name=SAMPLE_NAME,
            in_alignment_bam_chunk_files=alignment_chunk_bam_files,
            in_map_cores=MAP_CORES,
            in_map_disk=MAP_DISK,
            in_map_mem=MAP_MEM
    }
    
    if (REFERENCE_PREFIX != "") {
        # strip all the GRCh38's off our path list file.  we need them for surject as they are in the path
        # but fixBAMContigNaming above stripped them, so we don't need them downstream
        call fixPathNames {
            input:
                in_path_file=pipeline_path_list_file,
                in_prefix_to_strip=REFERENCE_PREFIX,
        }
    }
    File properly_named_path_list_file = select_first([fixPathNames.fixed_path_list_file, pipeline_path_list_file])
             
    # Split merged alignment by contigs list
    call splitBAMbyPath { 
        input:
            in_sample_name=SAMPLE_NAME,
            in_merged_bam_file=mergeAlignmentBAMChunks.merged_bam_file,
            in_merged_bam_file_index=mergeAlignmentBAMChunks.merged_bam_file_index,
            in_path_list_file=properly_named_path_list_file,
            in_map_cores=MAP_CORES,
            in_map_disk=MAP_DISK,
            in_map_mem=MAP_MEM
    }

    ##
    ## Prepare each BAM
    ##
    scatter (bam_and_index_for_path in zip(splitBAMbyPath.bam_contig_files, splitBAMbyPath.bam_contig_files_index)) {
        ## Evantually shift and realign reads
        if (LEFTALIGN_BAM){
            # Just left-shift each read individually
            call leftShiftBAMFile {
                input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=bam_and_index_for_path.left,
                in_reference_file=reference_file,
                in_reference_index_file=reference_index_file,
                in_call_disk=POSTPROCESS_DISK
            }
            # This tool can't make an index itself so we need to re-index the BAM
            call indexBAMFile {
                input:
                in_sample_name=SAMPLE_NAME,
                in_bam_file=leftShiftBAMFile.left_shifted_bam,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
            }
        }
        if (REALIGN_INDELS) {
            File forrealign_bam = select_first([leftShiftBAMFile.left_shifted_bam, bam_and_index_for_path.left])
            File forrealign_index = select_first([indexBAMFile.bam_index, bam_and_index_for_path.right])
            # Do indel realignment
            call runGATKRealignerTargetCreator {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_file=forrealign_bam,
                    in_bam_index_file=forrealign_index,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_reference_dict_file=reference_dict_file,
                    in_call_disk=POSTPROCESS_DISK
            }
            if (REALIGNMENT_EXPANSION_BASES != 0) {
                # We want the realignment targets to be wider
                call widenRealignmentTargets {
                    input:
                        in_target_bed_file=runGATKRealignerTargetCreator.realigner_target_bed,
                        in_reference_index_file=reference_index_file,
                        in_expansion_bases=REALIGNMENT_EXPANSION_BASES,
                        in_call_disk=POSTPROCESS_DISK
                }
            }
            File target_bed_file = select_first([widenRealignmentTargets.output_target_bed_file, runGATKRealignerTargetCreator.realigner_target_bed])
            call runAbraRealigner {
                input:
                    in_sample_name=SAMPLE_NAME,
                    in_bam_file=forrealign_bam,
                    in_bam_index_file=forrealign_index,
                    in_target_bed_file=target_bed_file,
                    in_reference_file=reference_file,
                    in_reference_index_file=reference_index_file,
                    in_call_disk=POSTPROCESS_DISK
            }
        }
        File processed_bam = select_first([runAbraRealigner.indel_realigned_bam, leftShiftBAMFile.left_shifted_bam, bam_and_index_for_path.left])
        File processed_bam_index = select_first([runAbraRealigner.indel_realigned_bam_index, indexBAMFile.bam_index, bam_and_index_for_path.right])
    }
    
    if (OUTPUT_SINGLE_BAM){
        call mergeAlignmentBAMChunks as mergeBAM {
            input:
                in_sample_name=SAMPLE_NAME,
                in_alignment_bam_chunk_files=processed_bam,
                in_map_cores=MAP_CORES,
                in_map_disk=MAP_DISK,
                in_map_mem=MAP_MEM
        }
    }
    
    if (OUTPUT_CALLING_BAMS){
        Array[File] calling_bams = processed_bam
        Array[File] calling_bam_indexes = processed_bam_index
    }
    
    if (OUTPUT_GAF){
        scatter (gam_chunk_file in runVGGIRAFFE.chunk_gam_file) {
            call convertGAMtoGAF {
                input:
                in_xg_file=XG_FILE,
                in_gam_file= gam_chunk_file,
                in_vg_container=VG_CONTAINER,
                in_cores=MAP_CORES,
                in_disk=MAP_DISK,
                in_mem=MAP_MEM
            }
        }
        call mergeGAF {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gaf_chunk_files=convertGAMtoGAF.output_gaf,
            in_vg_container=VG_CONTAINER,
            in_disk=2*MAP_DISK
        }
    }

    if (OUTPUT_GAM){
        call mergeGAM {
            input:
            in_sample_name=SAMPLE_NAME,
            in_gam_chunk_files=runVGGIRAFFE.chunk_gam_file,
            in_vg_container=VG_CONTAINER,
            in_disk=2*MAP_DISK
        }
    }
    
    output {
        File? output_bam = mergeBAM.merged_bam_file
        File? output_bam_index = mergeBAM.merged_bam_file_index
        File? output_gaf = mergeGAF.output_merged_gaf
        File? output_gam = mergeGAM.output_merged_gam
        Array[File]? output_calling_bams = calling_bams
        Array[File]? output_calling_bam_indexes = calling_bam_indexes
    }   
}

########################
### TASK DEFINITIONS ###
########################

task convertCRAMtoFASTQ {
    input {
	File? in_cram_file
        File? in_ref_file
        File? in_ref_index_file
	Int in_cores
    }
    Int half_cores = in_cores / 2
    Int disk_size = round(5 * size(in_cram_file, 'G')) + 50
    command <<<
    # Set the exit code of a pipeline to that of the rightmost command
    # to exit with a non-zero status, or zero if all commands of the pipeline exit
    set -o pipefail
    # cause a bash script to exit immediately when a command fails
    set -e
    # cause the bash shell to treat unset variables as an error and exit immediately
    set -u
    # echo each line of the script to stdout so we can see what is happening
    set -o xtrace
    #to turn off echo do 'set +o xtrace'
    
    samtools collate -@ ~{half_cores} --reference ~{in_ref_file} -Ouf ~{in_cram_file} | samtools fastq -@ ~{half_cores} -1 reads.R1.fastq.gz -2 reads.R2.fastq.gz -0 reads.o.fq.gz -s reads.s.fq.gz -c 1 -N -
    >>>
    output {
        File output_fastq_1_file = "reads.R1.fastq.gz"
        File output_fastq_2_file = "reads.R2.fastq.gz"
    }
    runtime {
        preemptible: 2
        cpu: in_cores
        memory: "50 GB"
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }    
}

task splitReads {
    input {
        File in_read_file
        String in_pair_id
        Int in_reads_per_chunk
        Int in_split_read_cores
        Int in_split_read_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        CHUNK_LINES=$(( ~{in_reads_per_chunk} * 4 ))
        gzip -cd ~{in_read_file} | split -l $CHUNK_LINES --filter='pigz -p ~{in_split_read_cores} > ${FILE}.fq.gz' - "fq_chunk_~{in_pair_id}.part."
    >>>
    output {
        Array[File] output_read_chunks = glob("fq_chunk_~{in_pair_id}.part.*")
    }
    runtime {
        preemptible: 2
        time: 120
        cpu: in_split_read_cores
        memory: "2 GB"
        disks: "local-disk " + in_split_read_disk + " SSD"
        docker: "quay.io/glennhickey/pigz:2.3.1"
    }
}

task extractSubsetPathNames {
    input {
        File in_xg_file
        String in_vg_container
        Int in_extract_disk
        Int in_extract_mem
    }

    command {
        set -eux -o pipefail

        vg paths \
            --list \
            --xg ${in_xg_file} > path_list.txt

        grep -v _decoy path_list.txt | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM > path_list.sub.txt
    }
    output {
        File output_path_list_file = "path_list.sub.txt"
    }
    runtime {
        preemptible: 2
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: in_vg_container
    }
}

task extractReference {
    input {
        File in_xg_file
        File in_path_list_file
        String in_vg_container
        Int in_extract_disk
        Int in_extract_mem
    }

    command {
        set -eux -o pipefail

        # Subset to just the paths we care about (may be the whole file) so we
        # get a good dict with just those paths later
        vg paths \
           --extract-fasta \
           -p ${in_path_list_file} \
           --xg ${in_xg_file} > ref.fa
    }
    output {
        File reference_file = "ref.fa"
    }
    runtime {
        preemptible: 2
        memory: in_extract_mem + " GB"
        disks: "local-disk " + in_extract_disk + " SSD"
        docker: in_vg_container
    }
}

task runVGGIRAFFE {
    input {
        File in_left_read_pair_chunk_file
        File in_right_read_pair_chunk_file
        File in_xg_file
        File in_gbwt_file
        File in_ggbwt_file
        File in_dist_file
        File in_min_file
        String in_vg_container
        String in_giraffe_options
        String in_sample_name
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        READ_CHUNK_ID=($(ls ~{in_left_read_pair_chunk_file} | awk -F'.' '{print $(NF-2)}'))
        vg giraffe \
          --progress \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --sample "~{in_sample_name}" \
          ~{in_giraffe_options} \
          --output-format gam \
          -f ~{in_left_read_pair_chunk_file} -f ~{in_right_read_pair_chunk_file} \
          -x ~{in_xg_file} \
          -H ~{in_gbwt_file} \
          -g ~{in_ggbwt_file} \
          -d ~{in_dist_file} \
          -m ~{in_min_file} \
          -t ~{in_map_cores} > ~{in_sample_name}.${READ_CHUNK_ID}.gam
    >>>
    output {
        File chunk_gam_file = glob("*gam")[0]
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }
}

task surjectGAMtoBAM {
    input {
        File in_gam_file
        File in_xg_file
        File in_path_list_file
        String in_sample_name
        Int in_max_fragment_length
        String in_vg_container
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }
    String out_prefix = basename(in_gam_file, ".gam") 
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        vg surject \
          -F ~{in_path_list_file} \
          -x ~{in_xg_file} \
          -t ~{in_map_cores} \
          --bam-output \
          --sample ~{in_sample_name} \
          --read-group "ID:1 LB:lib1 SM:~{in_sample_name} PL:illumina PU:unit1" \
          --prune-low-cplx \
          --interleaved --max-frag-len ~{in_max_fragment_length} \
          ~{in_gam_file} > ~{out_prefix}.bam
    >>>
    output {
        File chunk_bam_file = "~{out_prefix}.bam"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: in_vg_container
    }

}

task convertGAMtoGAF {
    input {
        File in_xg_file
        File in_gam_file
        String in_vg_container
        Int in_cores
        Int in_disk
        String in_mem
    }
    String out_prefix = basename(in_gam_file, ".gam") 
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        vg convert -G ~{in_gam_file} ~{in_xg_file} | gzip > ~{out_prefix}.gaf.gz
    >>>
    output {
        File output_gaf = "~{out_prefix}.gaf.gz"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: in_mem + " GB"
        cpu: in_cores
        disks: "local-disk " + in_disk + " SSD"
        docker: in_vg_container
    }
}

task mergeGAF {
    input {
        String in_sample_name
        Array[File] in_gaf_chunk_files
        String in_vg_container
        Int in_disk
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        cat ~{sep=" " in_gaf_chunk_files} > ~{in_sample_name}.gaf.gz
    >>>
    output {
        File output_merged_gaf = "~{in_sample_name}.gaf.gz"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: "6GB"
        cpu: 1
        disks: "local-disk " + in_disk + " SSD"
        docker: in_vg_container
    }
}

task mergeGAM {
    input {
        String in_sample_name
        Array[File] in_gam_chunk_files
        String in_vg_container
        Int in_disk
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        cat ~{sep=" " in_gam_chunk_files} > ~{in_sample_name}.gam
    >>>
    output {
        File output_merged_gam = "~{in_sample_name}.gam"
    }
    runtime {
        preemptible: 2
        time: 300
        memory: "6GB"
        cpu: 1
        disks: "local-disk " + in_disk + " SSD"
        docker: in_vg_container
    }
}

task sortBAMFile {
    input {
        String in_sample_name
        File in_bam_chunk_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        samtools sort \
          --threads ~{in_map_cores} \
          ~{in_bam_chunk_file} \
          -O BAM > ~{in_sample_name}.positionsorted.bam
    >>>
    output {
        File sorted_chunk_bam = "~{in_sample_name}.positionsorted.bam"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task indexBAMFile {
    input {
        String in_sample_name
        File in_bam_file
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        # Never use Samtools 1.4.1 here! See https://github.com/samtools/samtools/issues/687
        samtools index \
          ~{in_bam_file} \
          index.bai
    >>>
    output {
        File bam_index = "index.bai"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: in_map_mem + " GB"
        cpu: 1
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task fixBAMContigNaming {
    input {
        File in_bam_file
        File in_ref_dict
        String in_prefix_to_strip
        Int in_map_cores
        Int in_map_disk = 5 * round(size(in_bam_file, 'G')) + 20
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        # patch the SQ fields from the dict into a new header
        samtools view -H ~{in_bam_file} | grep ^@HD > new_header.sam
        grep ^@SQ ~{in_ref_dict} | awk '{print $1 "\t" $2 "\t" $3}' >> new_header.sam
        samtools view -H ~{in_bam_file}  | grep -v ^@HD | grep -v ^@SQ >> new_header.sam

        OUTBAM=($(basename ~{in_bam_file} | sed "s/~{in_prefix_to_strip}//g"))
        mkdir fixed
        # insert the new header, and strip all instances of the prefix
        cat <(cat new_header.sam) <(samtools view ~{in_bam_file}) | \
            sed -e "s/~{in_prefix_to_strip}//g" | \
            samtools view --threads ~{in_map_cores} -O BAM > fixed/$OUTBAM
        
        samtools index fixed/$OUTBAM
    >>>
    output {
        File fixed_bam_file = glob("fixed/*bam")[0]
        File fixed_bam_index_file = glob("fixed/*bai")[0]
    }
    runtime {
        preemptible: 2
        time: 90
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task fixPathNames {
    input {
        File in_path_file
        String in_prefix_to_strip
     }

     command <<<
        sed -e "s/~{in_prefix_to_strip}//g" ~{in_path_file}  > "fixed_names_file"
     >>>
    output {
        File fixed_path_list_file = "fixed_names_file"
    }
    runtime {
        preemptible: 2
        time: 90
        memory: 2 + " GB"
        cpu: 1
        disks: "local-disk " + 32 + " SSD"
        docker: "quay.io/cmarkello/samtools_picard@sha256:e484603c61e1753c349410f0901a7ba43a2e5eb1c6ce9a240b7f737bba661eb4"
    }
}

task mergeAlignmentBAMChunks {
    input {
        String in_sample_name
        Array[File] in_alignment_bam_chunk_files
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        samtools merge \
          -f -p -c --threads ~{in_map_cores} \
          ~{in_sample_name}_merged.positionsorted.bam \
          ~{sep=" " in_alignment_bam_chunk_files} \
        && samtools index \
          ~{in_sample_name}_merged.positionsorted.bam
    >>>
    output {
        File merged_bam_file = "~{in_sample_name}_merged.positionsorted.bam"
        File merged_bam_file_index = "~{in_sample_name}_merged.positionsorted.bam.bai"
    }
    runtime {
        preemptible: 2
        time: 240
        memory: 5 + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task splitBAMbyPath {
    input {
        String in_sample_name
        File in_merged_bam_file
        File in_merged_bam_file_index
        File in_path_list_file
        Int in_map_cores
        Int in_map_disk
        String in_map_mem
    }

    command <<<
        set -eux -o pipefail

        ln -s ~{in_merged_bam_file} input_bam_file.bam
        ln -s ~{in_merged_bam_file_index} input_bam_file.bam.bai

        while read -r contig; do
            samtools view \
              -@ ~{in_map_cores} \
              -h -O BAM \
              input_bam_file.bam ${contig} \
              -o ~{in_sample_name}.${contig}.bam \
            && samtools index \
              ~{in_sample_name}.${contig}.bam
        done < "~{in_path_list_file}"
    >>>
    output {
        Array[File] bam_contig_files = glob("~{in_sample_name}.*.bam")
        Array[File] bam_contig_files_index = glob("~{in_sample_name}.*.bam.bai")
    }
    runtime {
        preemptible: 2
        memory: in_map_mem + " GB"
        cpu: in_map_cores
        disks: "local-disk " + in_map_disk + " SSD"
        docker: "biocontainers/samtools@sha256:3ff48932a8c38322b0a33635957bc6372727014357b4224d420726da100f5470"
    }
}

task runGATKRealignerTargetCreator {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_reference_file
        File in_reference_index_file
        File in_reference_dict_file
        Int in_call_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command 
        # to exit with a non-zero status, or zero if all commands of the pipeline exit 
        set -o pipefail
        # cause a bash script to exit immediately when a command fails 
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately 
        set -u
        # echo each line of the script to stdout so we can see what is happening 
        set -o xtrace
        #to turn off echo do 'set +o xtrace' 

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s "~{in_reference_file}" reference.fa
        ln -f -s "~{in_reference_index_file}" reference.fa.fai
        # And the dict must be adjacent to both
        ln -f -s "~{in_reference_dict_file}" reference.dict

        java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
          --remove_program_records \
          -drf DuplicateRead \
          --disable_bam_indexing \
          -nt "32" \
          -R reference.fa \
          -L ${CONTIG_ID} \
          -I input_bam_file.bam \
          --out forIndelRealigner.intervals

        awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > ~{in_sample_name}.${CONTIG_ID}.intervals.bed
    >>>
    output {
        File realigner_target_bed = glob("*.bed")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 16
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "broadinstitute/gatk3@sha256:5ecb139965b86daa9aa85bc531937415d9e98fa8a6b331cb2b05168ac29bc76b"
    }
}

task widenRealignmentTargets {
    input {
        File in_target_bed_file
        File in_reference_index_file
        Int in_expansion_bases
        Int in_call_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'
        
        BASE_NAME=($(ls ~{in_target_bed_file} | rev | cut -f1 -d'/' | rev | sed s/.bed$//g))

        # Widen the BED regions, but don't escape the chromosomes
        bedtools slop -i "~{in_target_bed_file}" -g "~{in_reference_index_file}" -b "~{in_expansion_bases}" > "${BASE_NAME}.widened.bed"
    >>>
    output {
        File output_target_bed_file = glob("*.widened.bed")[0]
    }
    runtime {
        preemptible: 2
        memory: 4 + " GB"
        cpu: 1
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    }
}

task runAbraRealigner {
    input {
        String in_sample_name
        File in_bam_file
        File in_bam_index_file
        File in_target_bed_file
        File in_reference_file
        File in_reference_index_file
        Int in_call_disk
        Int memoryGb = 40
        Int threadCount = 16
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        ln -f -s ~{in_bam_file} input_bam_file.bam
        ln -f -s ~{in_bam_index_file} input_bam_file.bam.bai
        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai

        java -Xmx~{memoryGb}G -jar /opt/abra2/abra2.jar \
          --targets ~{in_target_bed_file} \
          --in input_bam_file.bam \
          --out ~{in_sample_name}.${CONTIG_ID}.indel_realigned.bam \
          --ref reference.fa \
          --index \
          --threads ~{threadCount}
    >>>
    output {
        File indel_realigned_bam = glob("~{in_sample_name}.*.indel_realigned.bam")[0]
        File indel_realigned_bam_index = glob("~{in_sample_name}.*.indel_realigned*bai")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: memoryGb + " GB"
        cpu: threadCount
        disks: "local-disk " + in_call_disk + " SSD"
        # This used to be docker: "dceoy/abra2:latest" but they moved the tag
        # and it stopped working. A known good version has been rehosted on
        # Quay in case Docker Hub deletes it.
        docker: "quay.io/adamnovak/dceoy-abra2@sha256:43d09d1c10220cfeab09e2763c2c5257884fa4457bcaa224f4e3796a28a24bba"
    }
}

task leftShiftBAMFile {
    input {
        String in_sample_name
        File in_bam_file
        File in_reference_file
        File in_reference_index_file
        Int in_call_disk
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        CONTIG_ID=($(ls ~{in_bam_file} | rev | cut -f1 -d'/' | rev | sed s/^~{in_sample_name}.//g | sed s/.bam$//g | sed s/.indel_realigned$//g | sed s/.left_shifted$//g))

        # Reference and its index must be adjacent and not at arbitrary paths
        # the runner gives.
        ln -f -s ~{in_reference_file} reference.fa
        ln -f -s ~{in_reference_index_file} reference.fa.fai
        
        bamleftalign \
            <~{in_bam_file} \
            >~{in_sample_name}.${CONTIG_ID}.left_shifted.bam \
            --fasta-reference reference.fa \
            --compressed
    >>>
    output {
        File left_shifted_bam = glob("~{in_sample_name}.*.left_shifted.bam")[0]
    }
    runtime {
        preemptible: 2
        time: 180
        memory: 20 + " GB"
        cpu: 1
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "biocontainers/freebayes:v1.2.0-2-deb_cv1"
    }
}

task concatClippedVCFChunks {
    input {
        String in_sample_name
        Array[File] in_clipped_vcf_chunk_files
        Int in_call_disk
        Int in_call_mem
    }

    command {
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        mkdir bcftools.tmp
        bcftools concat -n ${sep=" " in_clipped_vcf_chunk_files} | bcftools sort -T bcftools.tmp -O z -o ${in_sample_name}.vcf.gz - && bcftools index -t -o ${in_sample_name}.vcf.gz.tbi ${in_sample_name}.vcf.gz
    }
    output {
        File output_merged_vcf = "${in_sample_name}.vcf.gz"
        File output_merged_vcf_index = "${in_sample_name}.vcf.gz.tbi"
    }
    runtime {
        preemptible: 2
        time: 60
        memory: in_call_mem + " GB"
        disks: "local-disk " + in_call_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}
