vg\_wdl
---------------
Eric T Dawson, Mike Lin and Charles Markello, Jean Monlong, Adam Novak, Parsa Eskandar
MIT License, 2023

[Workflow Description Language (WDL)](https://software.broadinstitute.org/wdl/documentation/quickstart) scripts
for [vg](https://github.com/vgteam/vg) workflows.

- [Workflows](#Workflows)
- [Usage](#usage)
- [Testing locally](#Testing-locally)
- [Citation](#citation)
- [Contributing, Help, Bugs and Requests](#Contributing-Help-Bugs-and-Requests)

## Workflows

- **Giraffe-DeepVariant workflows**. Either the full Giraffe-DeepVariant workflow, or parts of it are available:
    - [Giraffe-DeepVariant workflow](#giraffe-deepvariant-workflow) to perform the full workflow: starting from
      FASTQs/CRAM, align reads to a pangenome and run [DeepVariant](https://github.com/google/deepvariant).
    - [Giraffe workflow](#giraffe-workflow) to map reads and produce BAMs ready to use
      by [DeepVariant](https://github.com/google/deepvariant).
    - [Giraffe-DeepVariant from GAF workflow](#giraffe-deepvariant-from-gaf-workflow) to project reads aligned to a
      pangenome (GAF), prepare them and run [DeepVariant](https://github.com/google/deepvariant).
- [Happy workflow](#happy-workflow) to evaluate small variants against a truthset
  using [hap.py](https://github.com/Illumina/hap.py)/[vcfeval](https://github.com/RealTimeGenomics/rtg-tools).
- [GAF to sorted GAM workflow](#gaf-to-sorted-gam-workflow) to convert a GAF into a sorted and indexed GAM. E.g. to use
  with the [sequenceTubeMap](https://github.com/vgteam/sequenceTubeMap).
- [Giraffe SV workflow](#Giraffe-SV-workflow) to map short reads to a pangenome and genotype SVs
  with [vg](https://github.com/vgteam/vg).
- [Haplotype Sampling workflow](#Haplotype-Sampling-workflow) to create a personalized pangenome using haplotype
  sampling
- [Map-call workflow](#Map-call-workflow) to map reads and call small variants [vg](https://github.com/vgteam/vg),
  DeepVariant and GATK (legacy?).
- [Map-call Pedigree workflow](#Map-call-Pedigree-workflow) to map reads and call variants in a pedigree
  with [vg](https://github.com/vgteam/vg) (legacy?).

See also the [Going further](#Going-further) section for more details on some aspects and HOW-TOs:

- [Path list](#Path-list)
- [Read realignment](#Read-realignment)
- [Using the HPRC pangenomes](#HPRC-pangenomes)
- [Reference prefix removal](#Reference-prefix-removal)
- [CRAM input](#CRAM-input)
- [Single-end reads](#Single-end-reads)
- [Interleaved reads](#Interleaved-reads)
- [Unmapped reads](#Unmapped-reads)
- [Reads chunking](#Reads-chunking)

### Giraffe-DeepVariant workflow

The full workflow to go from sequencing reads (FASTQs, CRAM) to small variant calls (VCF).

- workflow file: [workflows/giraffe_and_deepvariant.wdl](workflows/giraffe_and_deepvariant.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/GiraffeDeepVariant:gbz?tab=info)
- If you use this workflow, please cite [the HPRC preprint](#cite-HPRC).

Parameters (semi-auto-generated from the parameter_meta section):

- *INPUT_READ_FILE_1*: Input sample 1st read pair fastq.gz  
- *INPUT_READ_FILE_2*: Input sample 2nd read pair fastq.gz  
- *INPUT_CRAM_FILE*: Input CRAM file  
- *CRAM_REF*: Genome fasta file associated with the CRAM file  
- *CRAM_REF_INDEX*: Index of the fasta file associated with the CRAM file  
- *GBZ_FILE*: Path to .gbz index file  
- *DIST_FILE*: Path to .dist index file  
- *MIN_FILE*: Path to .min index file  
- *ZIPCODES_FILE*: (OPTIONAL) For chaining-based alignment, path to .zipcodes index file
- *HAPL_FILE*: (OPTIONAL) Path to .hapl file used in haplotype sampling
- *SAMPLE_NAME*: The sample name  
- *OUTPUT_GAF*: Should a GAF file with the aligned reads be saved? Default is 'true'.  
- *OUTPUT_SINGLE_BAM*: Should a single merged BAM file be saved? If yes, unmapped reads will be included and 'calling bams' (one per contig) won't be outputted by default. Default is 'false'.  
- *OUTPUT_CALLING_BAMS*: Should individual contig BAMs used for calling be saved? Default is the opposite of OUTPUT_SINGLE_BAM.  
- *OUTPUT_UNMAPPED_BAM*: Should an unmapped reads BAM be saved? Default is false.
- *PAIRED_READS*: Are the reads paired? Default is 'true'.
- *INTERLEAVED_READS*: Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'.
- *READS_PER_CHUNK*: Number of reads contained in each mapping chunk. Default 20,000,000.  
- *CONTIGS*: (OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index.  
- *PATH_LIST_FILE*: (OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths.  
- *REFERENCE_PREFIX*: Remove this off the beginning of path names in surjected BAM (set to match prefix in PATH_LIST_FILE)  
- *REFERENCE_FILE*: (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required if the graph does not contain all bases of the reference.  
- *REFERENCE_INDEX_FILE*: (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.  
- *REFERENCE_DICT_FILE*: (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths.  
- *HAPLOID_CONTIGS*: (OPTIONAL) Names of contigs in the reference (without REFERENCE_PREFIX) that are haploid in this sample (often chrX and chrY). Not compatible with DeepVariant 1.5.  
- *PAR_REGIONS_BED_FILE*: (OPTIONAL) BED file with pseudo-autosomal regions. Not compatible with DeepVariant 1.5.  
- *PRUNE_LOW_COMPLEXITY*: Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realignment. Default is 'true'.  
- *LEFTALIGN_BAM*: Whether or not to left-align reads in the BAM. Default is 'true'.  
- *REALIGN_INDELS*: Whether or not to realign reads near indels. Default is 'true'.  
- *REALIGNMENT_EXPANSION_BASES*: Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160.  
- *MIN_MAPQ*: Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is the DeepVariant default for the model type.  
- *MAX_FRAGMENT_LENGTH*: Maximum distance at which to mark paired reads properly paired. Default is 3000.  
- *GIRAFFE_PRESET*: (OPTIONAL) Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10)  
- *GIRAFFE_OPTIONS*: (OPTIONAL) Extra command line options for Giraffe mapper  
- *TRUTH_VCF*: Path to .vcf.gz to compare against  
- *TRUTH_VCF_INDEX*: Path to Tabix index for TRUTH_VCF  
- *EVALUATION_REGIONS_BED*: BED to evaluate against TRUTH_VCF on, where false positives will be counted  
- *RESTRICT_REGIONS_BED*: BED to restrict comparison against TRUTH_VCF to  
- *TARGET_REGION*: Contig or region to restrict evaluation to  
- *RUN_STANDALONE_VCFEVAL*: Whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)  
- *DV_MODEL_TYPE*: Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA.  
- *DV_MODEL_META*: .meta file for a custom DeepVariant calling model  
- *DV_MODEL_INDEX*: .index file for a custom DeepVariant calling model  
- *DV_MODEL_DATA*: .data-00000-of-00001 file for a custom DeepVariant calling model  
- *DV_MODEL_FILES*: Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format  
- *DV_MODEL_VARIABLES_FILES*: Array of files that need to go in a 'variables' subdirectory for a DV model  
- *DV_KEEP_LEGACY_AC*: Should DV use the legacy allele counter behavior? Default is 'true'. Should be 'false' for HiFi.  
- *DV_NORM_READS*: Should DV normalize reads itself? Default is 'false'. Should be 'true' for HiFi.  
- *OTHER_MAKEEXAMPLES_ARG*: Additional arguments for the make_examples step of DeepVariant  
- *DV_IS_1_7_OR_NEWER*: Flag to use DeepVariant 1.7+ command line syntax and recommended flags. Must be true if providing a DV 1.7+ Docker image, and false if providing an older one.  
- *DV_NO_GPU_DOCKER*: Container image to use when running DeepVariant for steps that don't benefit from GPUs  
- *DV_GPU_DOCKER*: Container image to use when running DeepVariant for steps that benefit from GPUs  
- *SPLIT_READ_CORES*: Number of cores to use when splitting the reads into chunks. Default is 8.
- *SPLIT_READ_MEM*: Memory, in GB, to use when splitting the reads into chunks. Default is 50.
- *MAP_CORES*: Number of cores to use when mapping the reads. Default is 16.
- *MAP_MEM*: Memory, in GB, to use when mapping the reads. Default is 120.
- *HAPLOTYPE_SAMPLING*: Whether or not to use haplotype sampling before running giraffe. Default is 'true'.
- *INDEX_MINIMIZER_WEIGHTED*: Whether to use weighted minimizer indexing with haplotype sampling. (Default: true)
- *INDEX_MINIMIZER_MEM*: Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)
- *KMER_COUNTING_MEM*: Memory, in GB, to use when counting kmers. (Default: 64)
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)
- *BAM_PREPROCESS_MEM*: Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20.
- *REALIGN_MEM*: Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower.
- *CALL_CORES*: Number of cores to use when calling variants. Default is 8.  
- *CALL_MEM*: Memory, in GB, to use when calling variants. Default is 50.  
- *EVAL_MEM*: Memory, in GB, to use when evaluating variant calls. Default is 60.  
- *VG_DOCKER*: Container image to use when running vg  
- *VG_GIRAFFE_DOCKER*: Alternate container image to use when running vg giraffe mapping  
- *VG_SURJECT_DOCKER*: Alternate container image to use when running vg surject  


Related
topics: [read realignment](#Read-realignment), [reference prefix removal](#Reference-prefix-removal), [CRAM input](#CRAM-input), [reads chunking](#Reads-chunking), [path list](#Path-list), [single-end reads](#Single-end-reads), [interleaved reads](#Interleaved-reads), [unmapped reads](#Unmapped-reads), [HPRC pangenomes](#HPRC-pangenomes).

[Test locally](#testing-locally) with:

```sh
miniwdl run --as-me workflows/giraffe_and_deepvariant.wdl -i params/giraffe_and_deepvariant.json
miniwdl run --as-me workflows/giraffe_and_deepvariant.wdl -i params/giraffe_and_deepvariant_single_end.json
miniwdl run --as-me workflows/giraffe_and_deepvariant.wdl -i params/giraffe_and_deepvariant_cram.json
```

### Giraffe workflow

Core VG Giraffe mapping, usable for [DeepVariant](https://github.com/google/deepvariant).
Reads are mapped to a pangenome with [vg giraffe](https://github.com/vgteam/vg) and pre-processed (e.g. indel
realignment).

- workflow file: [workflows/giraffe.wdl](workflows/giraffe.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/Giraffe:gbz?tab=info)
- If you use this workflow, please cite [the HPRC preprint](#cite-HPRC).

Parameters (semi-auto-generated from the parameter_meta section):

- *INPUT_READ_FILE_1*: Input sample 1st read pair fastq.gz
- *INPUT_READ_FILE_2*: Input sample 2nd read pair fastq.gz
- *INPUT_CRAM_FILE*: Input CRAM file
- *CRAM_REF*: Genome fasta file associated with the CRAM file
- *CRAM_REF_INDEX*: Index of the fasta file associated with the CRAM file
- *GBZ_FILE*: Path to .gbz index file
- *DIST_FILE*: Path to .dist index file
- *MIN_FILE*: Path to .min index file
- *SAMPLE_NAME*: The sample name
- *OUTPUT_SINGLE_BAM*: Should a single merged BAM file be saved? If yes, unmapped reads will be inluded and 'calling
  bams' (one per contig) won't be outputed. Default is 'true'.
- *OUTPUT_CALLING_BAMS*: Should individual contig BAMs be saved? Default is 'false'.
- *OUTPUT_GAF*: Should a GAF file with the aligned reads be saved? Default is 'false'.
- *PAIRED_READS*: Are the reads paired? Default is 'true'.
- *INTERLEAVED_READS*: Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'.
- *READS_PER_CHUNK*: Number of reads contained in each mapping chunk. Default 20 000 000.
- *PATH_LIST_FILE*: (OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If
  neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths.
- *CONTIGS*: (OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index.
- *REFERENCE_PREFIX*: Remove this off the beginning of path names in surjected BAM (set to match prefix in
  PATH_LIST_FILE)
- *REFERENCE_FILE*: (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required
  if the graph does not contain all bases of the reference.
- *REFERENCE_INDEX_FILE*: (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.
- *REFERENCE_DICT_FILE*: (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if
  REFERENCE_INDEX_FILE i
- *LEFTALIGN_BAM*: Whether or not to left-align reads in the BAM. Default is 'true'.
- *REALIGN_INDELS*: Whether or not to realign reads near indels. Default is 'true'.
- *REALIGNMENT_EXPANSION_BASES*: Number of bases to expand indel realignment targets by on either side, to free up read
  tails in slippery regions. Default is 160.
- *MAX_FRAGMENT_LENGTH*: Maximum distance at which to mark paired reads properly paired. Default is 3000.
- *GIRAFFE_OPTIONS*: (OPTIONAL) extra command line options for Giraffe mapper
- *SPLIT_READ_CORES*: Number of cores to use when splitting the reads into chunks. Default is 8.
- *MAP_CORES*: Number of cores to use when mapping the reads. Default is 16.
- *MAP_MEM*: Memory, in GB, to use when mapping the reads. Default is 120.
- *BAM_PREPROCESS_MEM*: Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20.
- *REALIGN_MEM*: Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower.
- *HAPLOTYPE_SAMPLING*: Whether or not to use haplotype sampling before running giraffe. Default is 'true'
- *DIPLOID*:Whether or not to use diploid sampling while doing haplotype sampling. Has to use with Haplotype_sampling=true. Default is 'true'
- *HAPL_FILE*: (OPTIONAL) Path to .hapl file used in haplotype sampling
- *R_INDEX_FILE*: (OPTIONAL) Path to .ri file used in haplotype sampling
- *KFF_FILE*: (OPTIONAL) Path to .kff file used in haplotype sampling
- *HAPLOTYPE_NUMBER*: Number of generated synthetic haplotypes used in haplotype sampling. (Default: 4)
- *INDEX_MINIMIZER_MEM*: Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)
- *KMER_COUNTING_MEM*: Memory, in GB, to use when counting kmers. (Default: 64)
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)

Related
topics: [read realignment](#Read-realignment), [reference prefix removal](#Reference-prefix-removal), [CRAM input](#CRAM-input), [reads chunking](#Reads-chunking), [path list](#Path-list), [single-end reads](#Single-end-reads), [unmapped reads](#Unmapped-reads), [HPRC pangenomes](#HPRC-pangenomes), [Haplotype Sampling](#Haplotype-Sampling-workflow).

[Test locally](#testing-locally) with:

```sh
miniwdl run --as-me workflows/giraffe.wdl -i params/giraffe.json
miniwdl run --as-me workflows/giraffe.wdl -i params/giraffe.singleended.json
miniwdl run --as-me workflows/giraffe.wdl -i params/giraffe.singleended.cram.json
miniwdl run --as-me workflows/giraffe.wdl -i params/giraffe_and_haplotype_sampling.json
```

### Giraffe-DeepVariant from GAF workflow

Surject a GAF and prepare the BAMs (e.g. fix names, indel realign), and call small variants
with [DeepVariant](https://github.com/google/deepvariant).

- workflow file: [workflows/giraffe_and_deepvariant_fromGAF.wdl](workflows/giraffe_and_deepvariant_fromGAF.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/GiraffeDeepVariantFromGAF:gbz?tab=info)
- If you use this workflow, please cite [the HPRC preprint](#cite-HPRC).

Parameters (semi-auto-generated from the parameter_meta section):

- *INPUT_GAF*: Input gzipped GAF file
- *GBZ_FILE*: Path to .gbz index file
- *SAMPLE_NAME*: The sample name
- *OUTPUT_SINGLE_BAM*: Should a single merged BAM file be saved? If yes, unmapped reads will be inluded and 'calling
  bams' (one per contig) won't be outputed. Default is 'true'.
- *PAIRED_READS*: Are the reads paired? Default is 'true'.
- *PATH_LIST_FILE*: (OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If
  neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths.
- *CONTIGS*: (OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index.
- *REFERENCE_PREFIX*: Remove this off the beginning of path names in surjected BAM (set to match prefix in
  PATH_LIST_FILE)
- *REFERENCE_FILE*: (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required
  if the graph does not contain all bases of the reference.
- *REFERENCE_INDEX_FILE*: (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.
- *REFERENCE_DICT_FILE*: (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if
  REFERENCE_INDEX_FILE is set
- *LEFTALIGN_BAM*: Whether or not to left-align reads in the BAM. Default is 'true'.
- *REALIGN_INDELS*: Whether or not to realign reads near indels. Default is 'true'.
- *REALIGNMENT_EXPANSION_BASES*: Number of bases to expand indel realignment targets by on either side, to free up read
  tails in slippery regions. Default is 160.
- *MIN_MAPQ*: Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right
  than wrong. Default is 1
- *MAX_FRAGMENT_LENGTH*: Maximum distance at which to mark paired reads properly paired. Default is 3000.
- *DV_MODEL_META*: (OPTIONAL) .meta file for a custom DeepVariant calling model
- *DV_MODEL_INDEX*: (OPTIONAL) .index file for a custom DeepVariant calling model
- *DV_MODEL_DATA*: (OPTIONAL) .data-00000-of-00001 file for a custom DeepVariant calling model
- *DV_KEEP_LEGACY_AC*: Should DV use the legacy allele counter behavior? If unspecified this is done, unless the model is responsible for the setting.
- *DV_NORM_READS*: Should DV normalize reads itself? If unspecified this is not done, unless the model is responsible for the setting.
- *OTHER_MAKEEXAMPLES_ARG*: Additional arguments for the make_examples step of DeepVariant
- *VG_CORES*: Number of cores to use when projecting the reads. Default is 16.
- *VG_MEM*: Memory, in GB, to use when projecting the reads. Default is 120.
- *CALL_CORES*: Number of cores to use when calling variants. Default is 8.
- *CALL_MEM*: Memory, in GB, to use when calling variants. Default is 50.

Related
topics: [read realignment](#Read-realignment), [reference prefix removal](#Reference-prefix-removal), [path list](#Path-list), [single-end reads](#Single-end-reads), [interleaved reads](#Interleaved-reads), [unmapped reads](#Unmapped-reads), [HPRC pangenomes](#HPRC-pangenomes).

[Test locally](#testing-locally) with:

```sh
miniwdl run --as-me workflows/giraffe_and_deepvariant_fromGAF.wdl -i params/giraffe_and_deepvariant_gaf.json
miniwdl run --as-me workflows/giraffe_and_deepvariant_fromGAF.wdl -i params/giraffe_and_deepvariant_gaf_single_end.json
```

### Happy workflow

Evaluation of the small variant calls using [hap.py](https://github.com/Illumina/hap.py).

- workflow file: [workflows/happy_evaluation.wdl](workflows/happy_evaluation.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/HappyEvaluation:gbz?tab=info)

Parameters (semi-auto-generated from the parameter_meta section):

- *VCF*: bgzipped VCF with variant calls
- *VCF_INDEX*: (Optional) If specified, use this tabix index for the VCF instead of indexing it
- *TRUTH_VCF*: bgzipped VCF with truthset
- *TRUTH_VCF_INDEX*: (Optional) If specified, use this index for the truth VCF instead of indexing it
- *REFERENCE_FILE*: Use this FASTA reference.
- *REFERENCE_INDEX_FILE*: (Optional) If specified, use this .fai index instead of indexing the reference file.
- *EVALUATION_REGIONS_BED*: (Optional) BED to restrict comparison against TRUTH_VCF to
- *REFERENCE_PREFIX*: (Optional) Remove this off the beginning of sequence names in the VCF

[Test locally](#testing-locally) with:

```sh
miniwdl run --as-me workflows/happy_evaluation.wdl -i params/happy_evaluation.json
```

### GAF to sorted GAM workflow

Currently, only GAM file can be sorted and indexed, for example to extract and subgraph and visualize, or use with
the [sequenceTubeMap](https://github.com/vgteam/sequenceTubeMap).
This workflow converts reads aligned to a pangenome in a GAF file to a sorted and indexed GAM file.

- workflow file: [workflows/sort_graph_aligned_reads.wdl](workflows/sort_graph_aligned_reads.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/sortGraphAlignedReads:gbz?tab=info)

Parameters (semi-auto-generated from the parameter_meta section):

- *GAF_FILE*: GAF file to convert and sort.
- *GBZ_FILE*: the GBZ index of the graph
- *SAMPLE_NAME*: (Optional) a sample name

Related topics: [HPRC pangenomes](#HPRC-pangenomes).

[Test locally](#testing-locally) with:

```
miniwdl run --as-me workflows/sort_graph_aligned_reads.wdl -i params/sort_graph_aligned_reads.gaf.json
```

### Giraffe SV workflow

Workflow for mapping short reads and genotyping the structural variants in a pangenome.

- workflow file: [workflows/vg_map_call_sv.wdl](workflows/vg_map_call_sv.wdl)
- parameter file: [params/vg_map_call_sv_test.json](params/vg_map_call_sv_test.json)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/vg_map_call_sv:svpack?tab=info)
- If you use this workflow, please cite [the Giraffe-SV article](#cite-giraffe-sv).

### Haplotype Sampling workflow

Workflow for creating a personalized pangenome with [haplotype sampling](https://github.com/vgteam/vg/wiki/Haplotype-Sampling).

- [workflow file](https://github.com/vgteam/vg_wdl/blob/master/workflows/haplotype_sampling.wdl)
- [parameter file](https://github.com/vgteam/vg_wdl/blob/master/params/haplotype_sampling.json)

Parameters  (semi-auto-generated from the parameter_meta section):
- *GBZ_FILE*: Path to .gbz index file
- *INPUT_READ_FILE_FIRST*: Input sample 1st read pair fastq.gz
- *INPUT_READ_FILE_SECOND*: Input sample 2st read pair fastq.gz
- *HAPL_FILE*: Path to .hapl file
- *DIST_FILE*: Path to .dist file
- *R_INDEX_FILE*: Path to .ri file
- *KFF_FILE*: Path to .kff file
- *OUTPUT_NAME_PREFIX*: Name of the output file (Default: haplotype_sampled_graph)
- *KMER_LENGTH*: Size of kmer using for sampling (Up to 31) (Default: 29)
- *CORES*: Number of cores to use with commands. (Default: 16)
- *KMER_COUNTING_MEM*: Memory, in GB, to use when counting kmers. (Default: 64)
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)
- *INDEX_MINIMIZER_MEM*: Memory, in GB, to use when making the minimizer index. (Default: 320)
- *WINDOW_LENGTH*: Window length used for building the minimizer index. (Default: 11)
- *SUBCHAIN_LENGTH*: Target length (in bp) for subchains. (Default: 10000)
- *HAPLOTYPE_NUMBER*: Number of generated synthetic haplotypes. (Default: 4)
- *PRESENT_DISCOUNT*: Multiplicative factor for discounting scores for present kmers. (Default: 0.9)
- *HET_ADJUST*: Additive term for adjusting scores for heterozygous kmers. (Default: 0.05)
- *ABSENT_SCORE*: Score for absent kmers. (Default: 0.8)
- *INCLUDE_REFERENCE*: Include reference paths and generic paths from the full graph in the sampled graph. (Default:
true)
- *SET_REFERENCE*: Name of single reference to include in sampled graph. (Default: all references)
- *DIPLOID*: Activate diploid sampling. (Default: true)
- *INDEX_MINIMIZER_K*: K-mer size of minimizer index to produce for sampled graph. Should be 29 for short read mapping and 31 for long read mapping. (Default: 29)
- *INDEX_MINIMIZER_W*: Window size of minimizer index to produce for sampled graph. Should be 11 for short read mapping and 50 for long read mapping. (Default: 11)
- *INDEX_MINIMIZER_WEIGHTED*: Whether to produce a weighted minimizer index for the sampled graph. (Default: true)
- *VG_DOCKER*: Container image to use when running vg.

[Test locally](#testing-locally) with:

```sh
miniwdl run --as-me workflows/haplotype_sampling.wdl -i params/haplotype_sampling.json
```

### Map-call workflow

- workflow file: [workflows/vg_multi_map_call.wdl](workflows/vg_multi_map_call.wdl)
- parameter file: [params/vg_multi_map_call.inputs_tiny.json](params/vg_multi_map_call.inputs_tiny.json)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/vg-pipeline-workingexample:master?tab=info)

### Map-call Pedigree workflow

- workflow file: [workflows/vg_trio_multi_map_call.wdl](workflows/vg_trio_multi_map_call.wdl)
- parameter file: [params/vg_trio_multi_map_call.inputs_tiny.json](params/vg_trio_multi_map_call.inputs_tiny.json)
- If you use this workflow, please cite the [Pedigree-VG article](#Cite-Pedigree-VG).

### Going further

See below more information
about: [read realignment](#Read-realignment), [reference prefix removal](#Reference-prefix-removal), [CRAM input](#CRAM-input), [reads chunking](#Reads-chunking), [path list](#Path-list), [single-end reads](#Single-end-reads), [interleaved reads](#Interleaved-reads), [unmapped reads](#Unmapped-reads), [HPRC pangenomes](#HPRC-pangenomes).

#### Read realignment

Once the reads are projected to a linear reference, we've noticed that realigning the reads can improve the variant
calling with [DeepVariant](https://github.com/google/deepvariant).
This helps mostly for the small insertions-deletions (indels).

The full realignment process involves:

1. Leftaligning the reads
   with [freebayes' `bamleftalign`](https://manpages.debian.org/testing/freebayes/bamleftalign.1.en.html).
    - Can be enabled/disabled with the `LEFTALIGN_BAM` parameter
2. Identifying regions to realign further
   with [GATK RealignerTargetCreator](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md).
3. Expand those regions with [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/slop.html).
    - Number of bases to expand controlled by the `REALIGNMENT_EXPANSION_BASES` parameter.
4. Realigning the reads in those regions with [ABRA2](https://github.com/mozack/abra2).

The last 3 steps can be enabled/disabled with the `REALIGN_INDELS` parameter.

Although it produces the best variant calls, these extra steps increase the computational resources (and cost) of the
workflow.
For a lighter run, switch off those two realignment steps and use [DeepVariant](https://github.com/google/deepvariant)'s
integrated realigner instead with:

- `LEFTALIGN_BAM=false`
- `REALIGN_INDELS=false`
- `DV_NORM_READS=true`

#### Reference prefix removal

The names of contigs/paths/haplotypes in pangenomes sometimes contains a prefix that we'd want to remove.
In the HPRC pangenomes, for example, the chromosomal contigs from GRCh38 are named *GRCh38.chr1*, etc.
In practice, we want to remove this prefix from the variant calls (VCFs), or reads aligned to that reference (BAMs).

This is controlled by the `REFERENCE_PREFIX` parameters in the workflows.
Setting `REFERENCE_PREFIX="GRCh38."` for example will ensure the VCFs/BAMs have *chr1*, etc. for contig names.

Because the pangenome uses them, **the prefix must still be present when specifying the paths to project the reads** too
though.
Hence, the `CONTIGS` and `PATH_LIST_FILE` must use the prefix.

However, **provided reference FASTAs or dictionary must not have the prefix**.
These could be FASTA or `.dict` files from the "official" reference genome or pre-computed for them, hence no prefix.
So, no prefix in `REFERENCE_FILE`, `REFERENCE_INDEX_FILE`, `REFERENCE_DICT_FILE`.

#### CRAM input

When the input is a CRAM file (`INPUT_CRAM_FILE`) instead of a pair of FASTQ
files (`INPUT_READ_FILE_1`/`INPUT_READ_FILE_2`), the user must also provide the appropriate reference FASTA to work with
that CRAM file with `CRAM_REF`, and its index with `CRAM_REF_INDEX`.

The CRAM file will be converted back to a pair of FASTQs, so it costs a little bit more to analyze CRAMs than FASTQs
currently.

#### Reads chunking

Sequencing reads are chunked to parallelize read mapping.
The amount of chunking is controlled by the `READS_PER_CHUNK` parameter which specify how many reads each chunk should
have.
For a WGS experiment, we use chunks of about 20M reads.

#### Path list

We might not always want to project the reads alignments to all the paths in the pangenome.
For example, we might only care about alignment to chromosomes and not alternate contigs.
Or there might be multiple sets of paths like in the CHM13-based HPRC pangenome which contains both reference paths for
CHM13 and GRCh38.
In that case, we can specify a list of paths to project the reads to using one of the following.

`PATH_LIST_FILE` is a file which lists the paths names, one per line.
For the HPRC pangenomes it looks like:

```txt
GRCh38.chr1
GRCh38.chr2
GRCh38.chr3
...etc
```

Otherwise, paths can be listed in the `CONTIGS` parameter as a list (WDL array).

#### Single-end reads

Workflows expect paired-end reads, but some workflows can also analyze single-end reads.

To use single-end reads:

- If providing FASTQs, only provide `INPUT_READ_FILE_1` (no `INPUT_READ_FILE_2`).
- Use `PAIRED_READS=false`

#### Interleaved reads

Some paired-end reads are stored in a single FASTQ file with the two reads of each pair interleaved.

To use interleaved paired-end reads:

- Only provide `INPUT_READ_FILE_1` (no `INPUT_READ_FILE_2`).
- Use `PAIRED_READS=true`
- Use `INTERLEAVED_READS=true`
- Ensure `READS_PER_CHUNK` is even

#### Unmapped reads

If including unmapped reads in the BAMs is important, make sure to switch on `OUTPUT_SINGLE_BAM=true` in
the [Giraffe-DeepVariant workflow](#Giraffe-DeepVariant-workflow) and [Giraffe workflow](#Giraffe-workflow).

#### HPRC pangenomes

We recommend using the filtered CHM13-based pangenome (freeze 1).
It contains both the CHM13 and GRCh38 reference paths.

Use the following indexes for the pangenome:

- `GBZ_FILE`: [hprc-v1.0-mc-chm13-minaf.0.1.gbz](https://storage.googleapis.com/hprc-pangenomes/hprc-v1.0-mc-chm13-minaf.0.1.gbz)
    - GBZ with the pangenome and haplotypes.
- `MIN_FILE`: [hprc-v1.0-mc-chm13-minaf.0.1.min](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/filtered/hprc-v1.0-mc-chm13-minaf.0.1.min)
    - Minimizer index.
- `DIST_FILE`: [hprc-v1.0-mc-chm13-minaf.0.1.dist](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/filtered/hprc-v1.0-mc-chm13-minaf.0.1.dist)
    - Distance index.

To project reads and call variants relative to the GRCh38 reference:

- `REFERENCE_PREFIX="GRCh38."`
- `PATH_LIST_FILE` containing *GRCh38.chr1*, *GRCh38.chr2*, etc. File available
  at [GRCh38.path_list.txt](https://storage.googleapis.com/hprc-pangenomes/GRCh38.path_list.txt)
- `REFERENCE_FILE`: [hg38.fa](https://storage.googleapis.com/hprc-pangenomes/hg38.fa)
- `REFERENCE_INDEX_FILE`: [hg38.fa.fai](https://storage.googleapis.com/hprc-pangenomes/hg38.fa.fai). Optional, the
  workflow will create it if necessary (for a small extra cost/time).
- `REFERENCE_DICT_FILE`: [hg38.dict](https://storage.googleapis.com/hprc-pangenomes/hg38.dict). Optional, the workflow
  will create it if necessary (for a small extra cost/time).

To project reads and call variants relative to the CHM13 reference:

- `REFERENCE_PREFIX="CHM13."`
- `PATH_LIST_FILE` containing *CHM13.chr1*, *CHM13.chr2*, etc. File available
  at [CHM13.path_list.txt](https://storage.googleapis.com/hprc-pangenomes/CHM13.path_list.txt)
- `REFERENCE_FILE`: [chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa](https://storage.googleapis.com/hprc-pangenomes/chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa)
- `REFERENCE_INDEX_FILE`: [chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa.fai](https://storage.googleapis.com/hprc-pangenomes/chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa.fai).
  Optional, the workflow will create it if necessary (for a small extra cost/time).
- `REFERENCE_DICT_FILE`: [chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.dict](https://storage.googleapis.com/hprc-pangenomes/chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.dict).
  Optional, the workflow will create it if necessary (for a small extra cost/time).

For earlier versions of DeepVariant (<1.5), models were retrained using reads aligned to the HPRC pangenomes.
The corresponding model files were deposited
at: [https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/PANGENOME_2022/DeepVariant/models/DEEPVARIANT_MC_Y1/](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/PANGENOME_2022/DeepVariant/models/DEEPVARIANT_MC_Y1/).
They can be passed to the workflows using the `DV_MODEL_META`, `DV_MODEL_INDEX`, and `DV_MODEL_DATA`.
Note that **it is not necessary to use custom models in the latest version of the workflows** as DeepVariant v1.5
includes default models suited for analyzing reads mapped to pangenomes (and projected back to a linear reference).

## Usage

### Dockstore

The workflows that were deposited on [Dockstore](https://dockstore.org/) can be launched
using [its command line](https://docs.dockstore.org/en/stable/launch-with/launch.html) or on platform
like [Terra](https://app.terra.bio/).

### Using miniwdl

[Install miniwdl](https://miniwdl.readthedocs.io/en/latest/getting_started.html#install-miniwdl), for example,
with `pip`:

```sh
pip3 install miniwdl
```

Clone this repo somewhere with `git clone https://github.com/vgteam/vg_wdl.git`

Run a workflow using:

```
miniwdl run /path/to/vg_wdl/workflows/WORKFLOW.wdl -i your-inputs.json
```

To modify the input parameters, edit the input `.json` with the necessary changes.

### Using Cromwell

[Cromwell](https://cromwell.readthedocs.io/en/stable/) can be run WDL workflows with:

```sh
java -jar $CROMWELL_JAR run workflow.wdl -i inputs.json
```

Where *CROMWELL_JAR* points at the Cromwell jar
downloaded [their release page](https://github.com/broadinstitute/cromwell/releases/), for example set
with `CROMWELL_JAR=/path/to/cromwell-<whatever>.jar` in your shell.

To run one of the workflows in this repo, clone the repo somewhere with `git clone https://github.com/vgteam/vg_wdl.git`
and run the desired workflow `.wdl` file:

```sh
java -jar $CROMWELL_JAR run /path/to/vg_wdl/workflows/WORKFLOW.wdl -i inputs.json
```

### Docker Containers

WDL needs the runtime Docker image to be present online (e.g. Dockerhub).
[Cromwell](#using-cromwell)/[miniwdl](#using-miniwdl) will pull those images automatically.
VG images are available at [quay.io](https://quay.io/repository/vgteam/vg?tab=tags) and can be pulled with:

```
docker pull quay.io/vgteam/vg:v1.44.0
```

Specific versions can be specified like above for version `v1.44.0`.

## Testing locally

To test the workflow locally, e.g. on the [small simulated dataset](tests/small_sim_graph), you can run it with Cromwell
or miniwdl (see [Usage](#usage)).
So, from the root of this repo, run something like:

```sh
java -jar $CROMWELL_JAR run workflows/WORKFLOW.wdl -i params/INPUTS.json
## or
miniwdl run --as-me workflows/WORKFLOW.wdl -i params/INPUTS.json
```

[Miniwdl](#using-miniwdl) might be slightly more useful when developing/testing a WDL because is catches errors in WDL
syntax faster, and is a bit more explicit about them.

## Citation

### Cite HPRC

If you use the Giraffe-DeepVariant workflows, please cite
the [HPRC preprint](https://www.biorxiv.org/content/10.1101/2022.07.09.499321v1):

```
Liao, Asri, Ebler, et al. A Draft Human Pangenome Reference. preprint, bioRxiv 2022; doi: https://doi.org/10.1101/2022.07.09.499321 
```

### Cite Giraffe-SV

If you use the SV genotyping workflow with vg giraffe, please
cite [this article](https://doi.org/10.1126/science.abg8871):

```
Sirén, Monlong, Chang, Novak, Eizenga, et al. Pangenomics Enables Genotyping of Known Structural Variants in 5202 Diverse Genomes. Science, vol. 374, no. 6574, Dec. 2021; doi: https://doi.org/10.1126/science.abg8871.
```

### Cite Pedigree-VG

If you use the pedigree-based workflow for rare variant discovery, please
cite [this article](https://pubmed.ncbi.nlm.nih.gov/35483961/):

```
Markello et al. A Complete Pedigree-Based Graph Workflow for Rare Candidate Variant Analysis. Genome Research, Apr. 2022; doi: https://doi.org/10.1101/gr.276387.121.
```

## Contributing, Help, Bugs and Requests

Please open an Issue on [GitHub](https://github.com/vgteam/vg_wdl/issues) for help, bug reports, or feature requests.
When doing so, please remember that vg\_wdl is open-source software made by a community of developers.
Please be considerate and support a positive environment.
