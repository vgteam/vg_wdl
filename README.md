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
    - [DeepVariant workflow](#deepvariant-workflow) to prepare already-mapped reads (BAM) and
      run [DeepVariant](https://github.com/google/deepvariant) on them.
- [Happy workflow](#happy-workflow) to evaluate small variants against a truthset
  using [hap.py](https://github.com/Illumina/hap.py)/[vcfeval](https://github.com/RealTimeGenomics/rtg-tools).
- [Aardvark workflow](#aardvark-workflow) to evaluate small variants against a truthset
  using [Aardvark](https://github.com/PacificBiosciences/aardvark).
- [Giraffe acceptance test workflow](#giraffe-acceptance-test-workflow) to compare two vg versions, by mapping and
  calling the same sample with each and evaluating both call sets against the same truth set.
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

The workflows above call [internal subworkflows](#internal-subworkflows) for the steps they share. Those are not meant
to be run on their own.

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
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/GiraffeDeepVariant:master?tab=info)
- If you use this workflow, please cite [the HPRC preprint](#cite-HPRC).

Parameters:

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
- *OUTPUT_SINGLE_BAM*: Should a single merged BAM file of reads used for calling be saved? If yes, unmapped reads will be included and 'calling bams' (one per contig) won't be outputted by default. Default is 'false'.
- *OUTPUT_CALLING_BAMS*: Should individual contig BAMs used for calling be saved? Default is the opposite of OUTPUT_SINGLE_BAM.
- *OUTPUT_UNMAPPED_BAM*: Should an unmapped reads BAM be saved? Default is false.
- *PAIRED_READS*: Are the reads paired? Default is 'true'.
- *INTERLEAVED_READS*: Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'.
- *READS_PER_CHUNK*: Number of reads contained in each mapping chunk. Default 20 million.
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
- *EVALUATION_REGIONS_BED*: BED to evaluate against TRUTH_VCF on, where false positives will be counted. Required when EVALUATE_WITH_AARDVARK is set.
- *EVALUATE_WITH_AARDVARK*: Should the calls be compared to TRUTH_VCF with Aardvark instead of hap.py? Default is 'false'.
- *STRATIFICATION_ARCHIVE*: (OPTIONAL) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files) to break the results down by. Only used when EVALUATE_WITH_AARDVARK is set.
- *RESTRICT_REGIONS_BED*: BED to restrict comparison against TRUTH_VCF to
- *TARGET_REGION*: Contig or region to restrict evaluation to
- *RUN_STANDALONE_VCFEVAL*: Whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)
- *DV_MODEL_TYPE*: Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA.
- *DV_MODEL_META*: .meta file for a custom DeepVariant calling model
- *DV_MODEL_INDEX*: .index file for a custom DeepVariant calling model
- *DV_MODEL_DATA*: .data-00000-of-00001 file for a custom DeepVariant calling model
- *DV_MODEL_FILES*: Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format
- *DV_MODEL_VARIABLES_FILES*: Array of files that need to go in a 'variables' subdirectory for a DV model
- *DV_PANGENOME_GBZ*: (OPTIONAL) Path to a pangenome graph in GBZ format for pangenome-aware DV.
- *DV_PANGENOME_IMAGE_HEIGHT*: (OPTIONAL) Height of the pangenome part of the pileup images for pangenome-aware DV. It will be used only if DV_PANGENOME_GBZ is set. If DV_PANGENOME_HAPLOTYPE_SAMPLING is done by this workflow and this is not set, it defaults to DV_PANGENOME_HAPLOTYPE_NUMBER + 5, DeepVariant's convention for a graph with that many haplotypes. If passing in an already-sampled DV_PANGENOME_GBZ instead, set this explicitly to (haplotype count + 5); leaving it unset then gets DeepVariant's own default, which is tuned for the un-sampled reference pangenome.
- *DV_PANGENOME_SHARED_MEMORY_SIZE_GB*: (OPTIONAL) Size of the shared memory segment in GB for loading pangenome in DeepVariant. It will be used only if PANGENOME_GBZ is set.
- *DV_PANGENOME_REFERENCE_PREFIX*: (OPTIONAL) Prefix on chromosome names in the pangenome GBZ (like 'GRCh38.') that isn't on the corresponding names in the BAM, analogous to REFERENCE_PREFIX but for the pangenome reference instead of the calling reference. Empty by default.
- *DV_PANGENOME_REF_NAME*: (OPTIONAL) The name of the reference to keep in the pangenome gbz file for pangenome-aware DV; all other reference-sense paths are removed before calling. Required if DV_PANGENOME_GBZ is set.
- *DV_PANGENOME_HAPLOTYPE_SAMPLING*: Should haplotype sampling of DV_PANGENOME_GBZ be done before pangenome-aware DV calling? This is a separate round of sampling from HAPLOTYPE_SAMPLING, which (if used) samples GBZ_FILE before mapping. Default is 'false'.
- *DV_PANGENOME_DIPLOID_SAMPLING*: Should the DV_PANGENOME_HAPLOTYPE_SAMPLING round of haplotype sampling be done in diploid mode? Default is 'false'.
- *DV_PANGENOME_HAPLOTYPE_NUMBER*: Number of haplotypes to sample for DV_PANGENOME_HAPLOTYPE_SAMPLING. Also used, if DV_PANGENOME_IMAGE_HEIGHT is not set, to size the pangenome-aware DV pileup images, so set it to the actual haplotype count even when passing in an already-sampled DV_PANGENOME_GBZ. Default is 32.
- *DV_KEEP_LEGACY_AC*: Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads.
- *DV_NORM_READS*: Should DV normalize reads itself? If unspecified this is not done, unless set in the model.
- *OTHER_MAKEEXAMPLES_ARG*: Additional arguments for the make_examples step of DeepVariant
- *DV_USE_GPUS*: Should DeepVariant use GPUs for calling variants? Default is 'true'.
- *DV_NO_GPU_DOCKER*: Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+.
- *DV_GPU_DOCKER*: Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+.
- *SPLIT_READ_CORES*: Number of cores to use when splitting the reads into chunks. Default is 8.
- *SPLIT_READ_MEM*: Memory, in GB, to use when splitting the reads into chunks. Default is 50.
- *MAP_CORES*: Number of cores to use when mapping the reads. Default is 16.
- *MAP_MEM*: Memory, in GB, to use when mapping the reads. Default is 120.
- *HAPLOTYPE_SAMPLING*: Whether or not to use haplotype sampling before running giraffe. Default is 'true'.
- *INDEX_MINIMIZER_WEIGHTED*: Whether to use weighted minimizer indexing with haplotype sampling. (Default: true)
- *INDEX_MINIMIZER_MEM*: Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)
- *KMER_COUNTING_MEM*: Memory, in GB, to use when counting kmers. (Default: 64)
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 200)
- *BAM_PREPROCESS_MEM*: Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20.
- *REALIGN_MEM*: Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower.
- *CALL_CORES*: Number of cores to use when calling variants. Default is 8.
- *CALL_MEM*: Memory, in GB, to use when calling variants. Default is 50.
- *MAKE_EXAMPLES_CORES*: Number of cores to use when making DeepVariant examples. Default is CALL_CORES.
- *MAKE_EXAMPLES_MEM*: Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM.
- *EVAL_CORES*: Number of cores to use when evaluating variant calls. Default is 8.
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

[params/giraffe_and_deepvariant_pangenome.json](params/giraffe_and_deepvariant_pangenome.json) shows how to wire up
pangenome-aware DeepVariant; running it needs a
[pangenome-aware DeepVariant container](https://www.biorxiv.org/content/10.1101/2025.06.05.657102v1) instead of the
usual DeepVariant one.

### Giraffe workflow

Core VG Giraffe mapping, usable for [DeepVariant](https://github.com/google/deepvariant).
Reads are mapped to a pangenome with [vg giraffe](https://github.com/vgteam/vg) and pre-processed (e.g. indel
realignment).

- workflow file: [workflows/giraffe.wdl](workflows/giraffe.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/Giraffe:master?tab=info)
- If you use this workflow, please cite [the HPRC preprint](#cite-HPRC).

Parameters:

- *INPUT_READ_FILE_1*: Input sample 1st read pair fastq.gz or fastq
- *INPUT_READ_FILE_2*: Input sample 2nd read pair fastq.gz or fastq
- *INPUT_CRAM_FILE*: Input CRAM file to realign
- *CRAM_REF*: Genome fasta file associated with the CRAM file
- *CRAM_REF_INDEX*: Index of the fasta file associated with the CRAM file
- *INPUT_BAM_FILE*: Input BAM file to realign
- *READ_CHUNKS_1*: (OPTIONAL) Input reads to map (either all reads or read 1), already split. When used, INPUT_READ_FILE_1 is still used for haplotype sampling.
- *READ_CHUNKS_2*: (OPTIONAL) Input reads to map (read 2), in the same order as READ_CHUNKS_1. Only used with READ_CHUNKS_1, when the reads are paired and not interleaved.
- *GBZ_FILE*: Path to .gbz index file
- *DIST_FILE*: (OPTIONAL) Path to .dist index file for the graph that will actually be mapped against (the
  haplotype-sampled graph, if HAPLOTYPE_SAMPLING is used). Generated from that graph if not given.
- *MIN_FILE*: (OPTIONAL) Path to .min index file for the graph that will actually be mapped against (the
  haplotype-sampled graph, if HAPLOTYPE_SAMPLING is used). Generated from that graph if not given.
- *ZIPCODES_FILE*: (OPTIONAL) For chaining-based alignment, path to .zipcodes index file
- *SAMPLE_NAME*: The sample name
- *OUTPUT_SINGLE_BAM*: Should a single merged BAM file be saved? Default is 'true'.
- *OUTPUT_CALLING_BAMS*: Should individual contig BAMs be saved? Default is 'false'.
- *OUTPUT_GAF*: Should a GAF file with the aligned reads be saved? Default is 'false'.
- *OUTPUT_GAF_CHUNKS*: Should the unmerged GAF chunks be saved? Default is 'false'.
- *PAIRED_READS*: Are the reads paired? Default is 'true'.
- *INTERLEAVED_READS*: Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'.
- *READS_PER_CHUNK*: Number of reads contained in each mapping chunk. Default 20 million.
- *PATH_LIST_FILE*: (OPTIONAL) Text file where each line is a path name in the GBZ index, to use instead of CONTIGS. If
  neither is given, paths are extracted from the GBZ and subset to chromosome-looking paths. If using REFERENCE_PREFIX,
  contig names in here should have the prefix.
- *CONTIGS*: (OPTIONAL) Desired reference genome contigs, which are all paths in the GBZ index. If using
  REFERENCE_PREFIX, contig names in here should have the prefix.
- *REFERENCE_PREFIX*: Remove this off the beginning of path names in surjected BAM (set to match prefix in
  PATH_LIST_FILE)
- *REFERENCE_FILE*: (OPTIONAL) If specified, use this FASTA reference instead of extracting it from the graph. Required
  if the graph does not contain all bases of the reference. If using REFERENCE_PREFIX, contig names in here should not
  have the prefix.
- *REFERENCE_INDEX_FILE*: (OPTIONAL) If specified, use this .fai index instead of indexing the reference file. If using
  REFERENCE_PREFIX, contig names in here should not have the prefix.
- *REFERENCE_DICT_FILE*: (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths. Required if
  REFERENCE_INDEX_FILE is set. If using REFERENCE_PREFIX, contig names in here should not have the prefix. This is used
  in BAM processing and not for choosing contigs for the surjection, which uses PATH_LIST_FILE.
- *PRUNE_LOW_COMPLEXITY*: Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realignment. Default is 'true'.
- *LEFTALIGN_BAM*: Whether or not to left-align reads in the BAM. Default is 'true'.
- *REALIGN_INDELS*: Whether or not to realign reads near indels. Default is 'true'.
- *REALIGNMENT_EXPANSION_BASES*: Number of bases to expand indel realignment targets by on either side, to free up read
  tails in slippery regions. Default is 160.
- *MAX_FRAGMENT_LENGTH*: Maximum distance at which to mark paired reads properly paired. Default is 3000.
- *GIRAFFE_PRESET*: (OPTIONAL) Name of Giraffe mapper parameter preset to use (default, fast, hifi, or r10)
- *GIRAFFE_OPTIONS*: (OPTIONAL) extra command line options for Giraffe mapper
- *SPLIT_READ_CORES*: Number of cores to use when splitting the reads into chunks. Default is 8.
- *SPLIT_READ_MEM*: Memory, in GB, to use when splitting the reads into chunks. Default is 50.
- *MAP_CORES*: Number of cores to use when mapping the reads. Default is 16.
- *MAP_MEM*: Memory, in GB, to use when mapping the reads. Default is 120.
- *BAM_PREPROCESS_MEM*: Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20.
- *REALIGN_MEM*: Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower.
- *HAPLOTYPE_SAMPLING*: Whether or not to use haplotype sampling before running giraffe. Default is 'true'
- *DIPLOID*:Whether or not to use diploid sampling while doing haplotype sampling. Has to use with Haplotype_sampling=true. Default is 'true'
- *SET_REFERENCE*: (OPTIONAL) Name of the single reference to keep for haplotype sampling.
- *HAPL_FILE*: (OPTIONAL) Path to .hapl file used in haplotype sampling
- *R_INDEX_FILE*: (OPTIONAL) Path to .ri file used in haplotype sampling
- *KFF_FILE*: (OPTIONAL) Path to .kff file used in haplotype sampling
- *HAPLOTYPE_NUMBER*: Number of generated synthetic haplotypes used in haplotype sampling. (Default: 32)
- *INDEX_MINIMIZER_WEIGHTED*: Whether to use weighted minimizer indexing with haplotype sampling. (Default: true)
- *INDEX_MINIMIZER_MEM*: Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)
- *KMER_COUNTING_MEM*: Memory, in GB, to use when counting kmers. (Default: 64)
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)
- *OUTPUT_HAPL*: Whether or not to output the haplotype index (.hapl) created before haplotype sampling. This is useful if the same pangenome will be used for DeepVariant calling after mapping. Default is 'false'.
- *VG_DOCKER*: Container image to use when running vg
- *VG_GIRAFFE_DOCKER*: Alternate container image to use when running vg giraffe mapping
- *VG_SURJECT_DOCKER*: Alternate container image to use when running vg surject

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
with [DeepVariant](https://github.com/google/deepvariant). Given a truth set, the calls are also compared to it,
with [hap.py](https://github.com/Illumina/hap.py) or, if *EVALUATE_WITH_AARDVARK* is set,
with [Aardvark](https://github.com/PacificBiosciences/aardvark).

The GAF can be given whole, in which case it is split up to surject in parallel, or as chunks that are already split.

- workflow file: [workflows/giraffe_and_deepvariant_fromGAF.wdl](workflows/giraffe_and_deepvariant_fromGAF.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/GiraffeDeepVariantFromGAF:master?tab=info)
- If you use this workflow, please cite [the HPRC preprint](#cite-HPRC).

Parameters:

- *INPUT_GAF*: (OPTIONAL) Input gzipped GAF file, which is split up to surject in parallel. Give this or GAF_CHUNKS.
- *GAF_CHUNKS*: (OPTIONAL) Input gzipped GAF, already split into chunks that can be surjected in parallel. Give this or INPUT_GAF.
- *READS_PER_CHUNK*: Number of reads to put in each chunk when splitting INPUT_GAF. Unused if GAF_CHUNKS is given. Default 20 million.
- *GBZ_FILE*: Path to .gbz index file. Has to be the graph the reads were mapped to, since the alignments name its nodes.
- *SAMPLE_NAME*: The sample name
- *OUTPUT_SINGLE_BAM*: Should a single merged BAM file of reads used for calling be saved? If yes, unmapped reads will be included and 'calling
  bams' (one per contig) won't be outputted. Default is 'true'.
- *OUTPUT_CALLING_BAMS*: Should individual contig BAMs used for calling be saved? Default is the opposite of OUTPUT_SINGLE_BAM.
- *OUTPUT_UNMAPPED_BAM*: Should an unmapped reads BAM be saved? Default is false.
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
- *HAPLOID_CONTIGS*: (OPTIONAL) Names of contigs in the reference (without REFERENCE_PREFIX) that are haploid in this sample (often chrX and chrY). Not compatible with DeepVariant 1.5.
- *PAR_REGIONS_BED_FILE*: (OPTIONAL) BED file with pseudo-autosomal regions. Not compatible with DeepVariant 1.5.
- *PRUNE_LOW_COMPLEXITY*: Whether or not to remove low-complexity or short in-tail anchors when surjecting and force tail realignment. Default is 'true'.
- *LEFTALIGN_BAM*: Whether or not to left-align reads in the BAM. Default is 'true'.
- *REALIGN_INDELS*: Whether or not to realign reads near indels. Default is 'true'.
- *REALIGNMENT_EXPANSION_BASES*: Number of bases to expand indel realignment targets by on either side, to free up read
  tails in slippery regions. Default is 160.
- *MIN_MAPQ*: Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right
  than wrong. Default is 1
- *MAX_FRAGMENT_LENGTH*: Maximum distance at which to mark paired reads properly paired. Default is 3000.
- *SURJECT_OPTIONS*: Extra command line options for vg surject.
- *TRUTH_VCF*: (OPTIONAL) Path to .vcf.gz to compare the calls against. Evaluation only runs if this and TRUTH_VCF_INDEX are given.
- *TRUTH_VCF_INDEX*: (OPTIONAL) Path to Tabix index for TRUTH_VCF
- *EVALUATION_REGIONS_BED*: (OPTIONAL) BED to evaluate against TRUTH_VCF on, where false positives will be counted. Required when EVALUATE_WITH_AARDVARK is set.
- *EVALUATE_WITH_AARDVARK*: Should the calls be compared to TRUTH_VCF with Aardvark instead of hap.py? Default is 'false'.
- *STRATIFICATION_ARCHIVE*: (OPTIONAL) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files) to break the results down by. Only used when EVALUATE_WITH_AARDVARK is set.
- *RESTRICT_REGIONS_BED*: (OPTIONAL) BED to restrict comparison against TRUTH_VCF to
- *TARGET_REGION*: (OPTIONAL) Contig or region to restrict evaluation to
- *RUN_STANDALONE_VCFEVAL*: Whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)
- *DV_MODEL_TYPE*: Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA.
- *DV_MODEL_META*: (OPTIONAL) .meta file for a custom DeepVariant calling model
- *DV_MODEL_INDEX*: (OPTIONAL) .index file for a custom DeepVariant calling model
- *DV_MODEL_DATA*: (OPTIONAL) .data-00000-of-00001 file for a custom DeepVariant calling model
- *DV_MODEL_FILES*: Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format
- *DV_MODEL_VARIABLES_FILES*: Array of files that need to go in a 'variables' subdirectory for a DV model
- *PANGENOME_GBZ*: (OPTIONAL) Path to a pangenome graph in GBZ format for pangenome-aware DV. All reference-sense paths in it other than DV_PANGENOME_REF_NAME are removed before calling, since they are not part of this sample and would otherwise show up as uninformative extra tracks in the pileup images.
- *DV_PANGENOME_IMAGE_HEIGHT*: (OPTIONAL) Height of the pangenome part of the pileup images for pangenome-aware models. It will be used only if PANGENOME_GBZ is set. If DV_PANGENOME_HAPLOTYPE_SAMPLING is done by this workflow and this is not set, it defaults to DV_PANGENOME_HAPLOTYPE_NUMBER + 5, DeepVariant's convention for a graph with that many haplotypes. If passing in an already-sampled PANGENOME_GBZ instead, set this explicitly to (haplotype count + 5); leaving it unset then gets DeepVariant's own default, which is tuned for the un-sampled reference pangenome.
- *DV_PANGENOME_SHARED_MEMORY_SIZE_GB*: (OPTIONAL) Size of the shared memory segment in GB for loading pangenome in DeepVariant. It will be used only if PANGENOME_GBZ is set.
- *DV_PANGENOME_REFERENCE_PREFIX*: (OPTIONAL) Prefix on chromosome names in the pangenome GBZ (like 'GRCh38.') that isn't on the corresponding names in the BAM, analogous to REFERENCE_PREFIX but for the pangenome reference instead of the calling reference. Empty by default.
- *DV_PANGENOME_REF_NAME*: (OPTIONAL) The name of the reference to keep in the pangenome gbz file for pangenome-aware DV; all other reference-sense paths are removed before calling. Required if PANGENOME_GBZ is set.
- *DV_PANGENOME_HAPLOTYPE_SAMPLING*: Should haplotype sampling of PANGENOME_GBZ be done before pangenome-aware DV calling? Default is 'false'.
- *DV_PANGENOME_READS_FOR_SAMPLING_1*: (OPTIONAL) First input read file for haplotype sampling
- *DV_PANGENOME_READS_FOR_SAMPLING_2*: (OPTIONAL) Second input read file for haplotype sampling (if paired)
- *DV_PANGENOME_DIPLOID_SAMPLING*: Should haplotype sampling be done in diploid mode? Default is 'false'.
- *DV_PANGENOME_HAPLOTYPE_NUMBER*: Number of haplotypes to sample for haplotype sampling. Also used, if DV_PANGENOME_IMAGE_HEIGHT is not set, to size the pangenome-aware DV pileup images, so set it to the actual haplotype count even when passing in an already-sampled PANGENOME_GBZ. Default is 32.
- *DV_PANGENOME_HAPL_FILE*: (OPTIONAL) Path to .hapl file used in haplotype sampling
- *DV_PANGENOME_DIST_FILE*: (OPTIONAL) Path to .dist file used in haplotype sampling
- *DV_PANGENOME_HAPLOTYPE_SAMPLE_CORES*: Number of cores to use for haplotype sampling. Default is 16.
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling). (Default: 120)
- *DV_KEEP_LEGACY_AC*: Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads.
- *DV_NORM_READS*: Should DV normalize reads itself? If unspecified this is not done, unless set in the model.
- *DV_USE_GPUS*: Should DeepVariant use GPUs for calling variants? Default is 'true'.
- *DV_NO_GPU_DOCKER*: Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+.
- *DV_GPU_DOCKER*: Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+.
- *OTHER_MAKEEXAMPLES_ARG*: Additional arguments for the make_examples step of DeepVariant
- *VG_CORES*: Number of cores to use when projecting the reads. Default is 16.
- *VG_MEM*: Memory, in GB, to use when projecting the reads. Default is 120.
- *BAM_PREPROCESS_MEM*: Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20.
- *REALIGN_MEM*: Memory, in GB, to use for Abra indel realignment. Default is 40.
- *CALL_CORES*: Number of cores to use when calling variants. Default is 8.
- *CALL_MEM*: Memory, in GB, to use when calling variants. Default is 50.
- *MAKE_EXAMPLES_CORES*: Number of cores to use when making DeepVariant examples. Default is CALL_CORES.
- *MAKE_EXAMPLES_MEM*: Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM.
- *EVAL_CORES*: Number of cores to use when evaluating variant calls. Default is 8.
- *EVAL_MEM*: Memory, in GB, to use when evaluating variant calls. Default is 60.
- *VG_DOCKER*: Container image to use when running vg
- *VG_SURJECT_DOCKER*: (OPTIONAL) Alternate container image to use when running vg surject

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
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/HappyEvaluation:master?tab=info)

Parameters:

- *VCF*: bgzipped VCF with variant calls
- *VCF_INDEX*: (Optional) If specified, use this tabix index for the VCF instead of indexing it
- *TRUTH_VCF*: bgzipped VCF with truthset
- *TRUTH_VCF_INDEX*: (Optional) If specified, use this index for the truth VCF instead of indexing it
- *REFERENCE_FILE*: Use this FASTA reference.
- *REFERENCE_INDEX_FILE*: (Optional) If specified, use this .fai index instead of indexing the reference file.
- *EVALUATION_REGIONS_BED*: (Optional) BED to restrict comparison against TRUTH_VCF to
- *RESTRICT_REGIONS_BED*: BED to restrict comparison against TRUTH_VCF to
- *TARGET_REGION*: contig or region to restrict evaluation to
- *REFERENCE_PREFIX*: (Optional) Remove this off the beginning of sequence names in the VCF
- *REMOVE_HOM_REFS*: (Optional) Should homozygous ref calls be removed? (might help if hap.py segfaults). Default 'false'.
- *RUN_STANDALONE_VCFEVAL*: whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)
- *EVAL_CORES*: Number of cores to use when evaluating variant calls. Default is 8.
- *EVAL_MEM*: Memory, in GB, to use when evaluating variant calls. Default is 60.

[Test locally](#testing-locally) with:

```sh
miniwdl run --as-me workflows/happy_evaluation.wdl -i params/happy_evaluation.json
```

### DeepVariant workflow

Partial workflow to go from mapped reads (BAM) to small variant calls (VCF). Reads are pre-processed (e.g. indel
realignment). DeepVariant then calls small variants. This is the calling half of the
[Giraffe-DeepVariant workflow](#giraffe-deepvariant-workflow), for when the reads are already mapped.

Given a truth set, the calls are also compared to it. The comparison is done
with [hap.py](https://github.com/Illumina/hap.py) by default, or
with [Aardvark](https://github.com/PacificBiosciences/aardvark) if *EVALUATE_WITH_AARDVARK* is set. Aardvark needs to be
told where the truth set is complete, so *EVALUATION_REGIONS_BED* is required when using it.

- workflow file: [workflows/deepvariant.wdl](workflows/deepvariant.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/DeepVariant:master?tab=info)

Parameters:

- *MERGED_BAM_FILE*: The all-contigs sorted BAM to call with.
- *MERGED_BAM_FILE_INDEX*: The .bai index for the input BAM file
- *SAMPLE_NAME*: The sample name
- *OUTPUT_SINGLE_BAM*: Should a single merged BAM file of reads used for calling be saved? If yes, unmapped reads will be included and 'calling bams' (one per contig) won't be outputted by default. Default is 'false'.
- *OUTPUT_CALLING_BAMS*: Should individual contig BAMs used for calling be saved? Default is the opposite of OUTPUT_SINGLE_BAM.
- *OUTPUT_UNMAPPED_BAM*: Should an unmapped reads BAM be saved? Default is false.
- *CONTIGS*: Contig path names to use as PATH_LIST_FILE. Must be set if PATH_LIST_FILE is not.
- *PATH_LIST_FILE*: Text file where each line is a contig name to evaluate on. Must be set if CONTIGS is not.
- *REFERENCE_PREFIX*: Remove this off the beginning of path names to get contig names in the BAM (set to match prefix in PATH_LIST_FILE)
- *REFERENCE_PREFIX_ON_BAM*: If true, the REFERENCE_PREFIX is also on the sequence names in the BAM header and needs to be removed.
- *REFERENCE_FILE*: FASTA reference to call against.
- *REFERENCE_INDEX_FILE*: (OPTIONAL) If specified, use this .fai index instead of indexing the reference file.
- *REFERENCE_DICT_FILE*: (OPTIONAL) If specified, use this pre-computed .dict file of sequence lengths.
- *HAPLOID_CONTIGS*: (OPTIONAL) Names of contigs in the reference (without REFERENCE_PREFIX) that are haploid in this sample (often chrX and chrY). Not compatible with DeepVariant 1.5.
- *PAR_REGIONS_BED_FILE*: (OPTIONAL) BED file with pseudo-autosomal regions. Not compatible with DeepVariant 1.5.
- *LEFTALIGN_BAM*: Whether or not to left-align reads in the BAM. Default is 'true'. If true, all input reads, including secondaries, must have the read sequence given.
- *REALIGN_INDELS*: Whether or not to realign reads near indels. Default is 'true'. If true, all input reads must be in a read group.
- *REALIGNMENT_EXPANSION_BASES*: Number of bases to expand indel realignment targets by on either side, to free up read tails in slippery regions. Default is 160.
- *MIN_MAPQ*: Minimum MAPQ of reads to use for calling. 4 is the lowest at which a mapping is more likely to be right than wrong. Default is the DeepVariant default for the model type.
- *TRUTH_VCF*: Path to .vcf.gz to compare against
- *TRUTH_VCF_INDEX*: Path to Tabix index for TRUTH_VCF
- *EVALUATION_REGIONS_BED*: BED to evaluate against TRUTH_VCF on, where false positives will be counted. Required when EVALUATE_WITH_AARDVARK is set.
- *EVALUATE_WITH_AARDVARK*: Should the calls be compared to TRUTH_VCF with Aardvark instead of hap.py? Default is 'false'.
- *STRATIFICATION_ARCHIVE*: (OPTIONAL) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files) to break the results down by. Only used when EVALUATE_WITH_AARDVARK is set.
- *RESTRICT_REGIONS_BED*: BED to restrict comparison against TRUTH_VCF to
- *TARGET_REGION*: contig or region to restrict evaluation to
- *RUN_STANDALONE_VCFEVAL*: whether to run vcfeval on its own in addition to hap.py (can crash on some DeepVariant VCFs)
- *DV_MODEL_TYPE*: Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA.
- *DV_MODEL_META*: .meta file for a custom DeepVariant calling model
- *DV_MODEL_INDEX*: .index file for a custom DeepVariant calling model
- *DV_MODEL_DATA*: .data-00000-of-00001 file for a custom DeepVariant calling model
- *DV_MODEL_FILES*: Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format
- *DV_MODEL_VARIABLES_FILES*: Array of files that need to go in a 'variables' subdirectory for a DV model
- *PANGENOME_GBZ*: (OPTIONAL) Path to a pangenome graph in GBZ format for pangenome-aware DV. All reference-sense paths in it other than DV_PANGENOME_REF_NAME are removed before calling, since they are not part of this sample and would otherwise show up as uninformative extra tracks in the pileup images.
- *DV_PANGENOME_IMAGE_HEIGHT*: (OPTIONAL) Height of the pangenome part of the pileup images for pangenome-aware models. It will be used only if PANGENOME_GBZ is set. If DV_PANGENOME_HAPLOTYPE_SAMPLING is done by this workflow and this is not set, it defaults to DV_PANGENOME_HAPLOTYPE_NUMBER + 5, DeepVariant's convention for a graph with that many haplotypes. If passing in an already-sampled PANGENOME_GBZ instead, set this explicitly to (haplotype count + 5); leaving it unset then gets DeepVariant's own default, which is tuned for the un-sampled reference pangenome.
- *DV_PANGENOME_SHARED_MEMORY_SIZE_GB*: (OPTIONAL) Size of the shared memory segment in GB for loading pangenome in DeepVariant. It will be used only if PANGENOME_GBZ is set.
- *DV_PANGENOME_REFERENCE_PREFIX*: (OPTIONAL) Prefix on chromosome names in the pangenome GBZ (like 'GRCh38.') that isn't on the corresponding names in the BAM, analogous to REFERENCE_PREFIX but for the pangenome reference instead of the calling reference. Empty by default.
- *DV_PANGENOME_REF_NAME*: (OPTIONAL) The name of the reference to keep in the pangenome gbz file for pangenome-aware DV; all other reference-sense paths are removed before calling. Required if PANGENOME_GBZ is set.
- *DV_PANGENOME_HAPLOTYPE_SAMPLING*: Should haplotype sampling of PANGENOME_GBZ be done before pangenome-aware DV calling? Default is 'false'.
- *DV_PANGENOME_READS_FOR_SAMPLING_1*: (OPTIONAL) First input read file for haplotype sampling
- *DV_PANGENOME_READS_FOR_SAMPLING_2*: (OPTIONAL) Second input read file for haplotype sampling (if paired)
- *DV_PANGENOME_DIPLOID_SAMPLING*: Should haplotype sampling be done in diploid mode? Default is 'false'.
- *DV_PANGENOME_HAPLOTYPE_NUMBER*: Number of haplotypes to sample for haplotype sampling. Also used, if DV_PANGENOME_IMAGE_HEIGHT is not set, to size the pangenome-aware DV pileup images, so set it to the actual haplotype count even when passing in an already-sampled PANGENOME_GBZ. Default is 32.
- *DV_PANGENOME_HAPL_FILE*: (OPTIONAL) Path to .hapl file used in haplotype sampling
- *DV_PANGENOME_DIST_FILE*: (OPTIONAL) Path to .dist file used in haplotype sampling
- *DV_PANGENOME_R_INDEX_FILE*: (OPTIONAL) Path to .ri file used in haplotype sampling
- *DV_PANGENOME_KFF_FILE*: (OPTIONAL) Path to .kff file used in haplotype sampling
- *KMER_COUNTING_MEM*: Memory, in GB, to use when counting kmers. (Default: 64)
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling). (Default: 120)
- *DV_PANGENOME_HAPLOTYPE_SAMPLE_CORES*: Number of cores to use for haplotype sampling. Default is 16.
- *DV_KEEP_LEGACY_AC*: Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads.
- *DV_NORM_READS*: Should DV normalize reads itself? If unspecified this is not done, unless set in the model.
- *OTHER_MAKEEXAMPLES_ARG*: Additional arguments for the make_examples step of DeepVariant
- *DV_USE_GPUS*: Should DeepVariant use GPUs for calling variants? Default is 'true'.
- *DV_NO_GPU_DOCKER*: Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+.
- *DV_GPU_DOCKER*: Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+.
- *VG_DOCKER*: Container image to use when running vg. Only used for pangenome-aware DV's haplotype sampling and reference-removal steps.
- *BAM_PREPROCESS_MEM*: Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20.
- *REALIGN_MEM*: Memory, in GB, to use for Abra indel realignment. Default is 40.
- *CALL_CORES*: Number of cores to use when calling variants. Default is 8.
- *CALL_MEM*: Memory, in GB, to use when calling variants. Default is 50.
- *MAKE_EXAMPLES_CORES*: Number of cores to use when making DeepVariant examples. Default is CALL_CORES.
- *MAKE_EXAMPLES_MEM*: Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM.
- *EVAL_CORES*: Number of cores to use when evaluating variant calls. Default is 8.
- *EVAL_MEM*: Memory, in GB, to use when evaluating variant calls. Default is 60.

### Aardvark workflow

Evaluation of the small variant calls using [Aardvark](https://github.com/PacificBiosciences/aardvark).

- workflow file: [workflows/aardvark_evaluation.wdl](workflows/aardvark_evaluation.wdl)

Parameters:

- *QUERY_VCF*: bgzipped VCF with variant calls to evaluate
- *QUERY_VCF_INDEX*: (Optional) tabix index for QUERY_VCF; will be indexed if not provided
- *TRUTH_VCF*: bgzipped VCF with truthset
- *TRUTH_VCF_INDEX*: (Optional) tabix index for TRUTH_VCF; will be indexed if not provided
- *REFERENCE_FILE*: FASTA reference
- *REFERENCE_INDEX_FILE*: (Optional) .fai index; will be indexed if not provided
- *REGIONS_BED*: BED of regions to restrict comparison to
- *STRATIFICATION_ARCHIVE*: (Optional) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files)
- *SAMPLE_NAME*: Sample name, used to name the output directory/archive
- *THREADS*: Number of threads for aardvark compare. Default 16.
- *EVAL_MEM*: Memory, in GB, to use when evaluating variant calls. Default is 30.

[Test locally](#testing-locally) with:

```sh
miniwdl run --as-me workflows/aardvark_evaluation.wdl -i params/aardvark_evaluation.json
```

### Acceptance testing workflow

Workflow for comparing the calling accuracy for two different vg versions, a "candidate" and a "baseline".

Runs indexing, mapping, and surjection stages, followed by calling with [DeepVariant](https://github.com/google/deepvariant) and evaluation against a truth set with [Aardvark](https://github.com/PacificBiosciences/aardvark). You can use different vg versions or the same vg version for each stage; stages that don't need to run separately will be run once. By default, indexing is done once, and mapping and surjection are done independently.

Indexes can be provided to skip the indexing stage. If the candidate vg cannot use the baseline vg's indexes, set *CANDIDATE_SEPARATE_INDEXES*. The candidate run then uses the *CANDIDATE_\** index inputs, and anything not passed there is built with the candidate container. With haplotype sampling on, the sampled graph and its indexes count as indexes: they are made once with the baseline vg by default, and once per run when *CANDIDATE_SEPARATE_INDEXES* is set.

To compare two versions of `vg surject` while mapping with the baseline vg in both runs:

```json
{
  "AcceptanceTest.BASELINE_VG_DOCKER": "quay.io/vgteam/vg:v1.64.0",
  "AcceptanceTest.CANDIDATE_VG_DOCKER": "quay.io/vgteam/vg:CANDIDATE",
  "AcceptanceTest.CANDIDATE_VG_GIRAFFE_DOCKER": "quay.io/vgteam/vg:v1.64.0"
}
```

To compare two versions of `vg surject` while mapping with the *candidate* vg in both runs:

```json
{
  "AcceptanceTest.BASELINE_VG_DOCKER": "quay.io/vgteam/vg:v1.64.0",
  "AcceptanceTest.CANDIDATE_VG_DOCKER": "quay.io/vgteam/vg:CANDIDATE",
  "AcceptanceTest.BASELINE_VG_GIRAFFE_DOCKER": "quay.io/vgteam/vg:CANDIDATE",
}
```

To compare two versions of `vg giraffe` for mapping, while surjecting the same way in both runs, set *BASELINE_VG_SURJECT_DOCKER* to match the candidate, or *CANDIDATE_VG_SURJECT_DOCKER* to match the baseline.

Surjection settings can also be compared, using *BASELINE_VG_SURJECT_OPTIONS* and *CANDIDATE_VG_SURJECT_OPTIONS*.

Each run is evaluated against the truth set with Aardvark. The workflow produces the Aardvark summary, the Aardvark full output directory, and the VCF for each condition.

BAM and GAF read alignments can also be requested. If GAF alignments are shared between the two conditions, a unified GAF will be produced. Otherwise, two GAFS will be produced for the baseline and candidate conditions.

- workflow file: [workflows/acceptance_test.wdl](workflows/acceptance_test.wdl)

Parameters:

- *BASELINE_VG_DOCKER*: Container image to use when running vg for the baseline run, which is the known-good version to compare against
- *CANDIDATE_VG_DOCKER*: Container image to use when running vg for the candidate run, which is the version under test
- *BASELINE_VG_GIRAFFE_DOCKER*: (OPTIONAL) Container image to use when running vg giraffe mapping in the baseline run, instead of BASELINE_VG_DOCKER. If the same as the candidate Giraffe docker, mapping only runs once.
- *BASELINE_VG_SURJECT_DOCKER*: (OPTIONAL) Container image to use when running vg surject in the baseline run, instead of BASELINE_VG_DOCKER
- *CANDIDATE_VG_GIRAFFE_DOCKER*: (OPTIONAL) Container image to use when running vg giraffe mapping in the candidate run, instead of CANDIDATE_VG_DOCKER. If the same as the baseline Giraffe docker, mapping only runs once.
- *CANDIDATE_VG_SURJECT_DOCKER*: (OPTIONAL) Container image to use when running vg surject in the candidate run, instead of CANDIDATE_VG_DOCKER
- *BASELINE_VG_SURJECT_OPTIONS*: (OPTIONAL) Extra command line options for vg surject in the baseline run
- *CANDIDATE_VG_SURJECT_OPTIONS*: (OPTIONAL) Extra command line options for vg surject in the candidate run
- *INPUT_READ_FILE_1*: Input sample 1st read pair fastq.gz
- *INPUT_READ_FILE_2*: Input sample 2nd read pair fastq.gz
- *INPUT_CRAM_FILE*: Input CRAM file. Converted to FASTQ once and shared by both runs.
- *CRAM_REF*: Genome fasta file associated with the CRAM file
- *CRAM_REF_INDEX*: Index of the fasta file associated with the CRAM file
- *GBZ_FILE*: Path to .gbz index file. Used by both runs unless CANDIDATE_SEPARATE_INDEXES is set.
- *DIST_FILE*: (OPTIONAL) Path to .dist index file. Built with the baseline vg if not provided.
- *MIN_FILE*: (OPTIONAL) Path to .min index file. Built with the baseline vg if not provided.
- *ZIPCODES_FILE*: (OPTIONAL) For chaining-based alignment, path to .zipcodes index file matching MIN_FILE
- *HAPL_FILE*: (OPTIONAL) Path to .hapl file used in haplotype sampling
- *CANDIDATE_SEPARATE_INDEXES*: Should the candidate run get its own indexes instead of sharing the baseline run's? Set this when the two vg versions cannot use each other's indexes. Default is 'false'.
- *CANDIDATE_GBZ_FILE*: (OPTIONAL) Path to .gbz index file for the candidate run. Only used if CANDIDATE_SEPARATE_INDEXES is set; defaults to GBZ_FILE.
- *CANDIDATE_DIST_FILE*: (OPTIONAL) Path to .dist index file for the candidate run. Only used if CANDIDATE_SEPARATE_INDEXES is set; built with the candidate vg if not provided.
- *CANDIDATE_MIN_FILE*: (OPTIONAL) Path to .min index file for the candidate run. Only used if CANDIDATE_SEPARATE_INDEXES is set; built with the candidate vg if not provided.
- *CANDIDATE_ZIPCODES_FILE*: (OPTIONAL) Path to .zipcodes index file for the candidate run, matching CANDIDATE_MIN_FILE. Only used if CANDIDATE_SEPARATE_INDEXES is set.
- *CANDIDATE_HAPL_FILE*: (OPTIONAL) Path to .hapl file for the candidate run. Only used if CANDIDATE_SEPARATE_INDEXES is set.
- *SAMPLE_NAME*: The sample name
- *TRUTH_VCF*: Path to .vcf.gz of truth calls to evaluate both runs against
- *TRUTH_VCF_INDEX*: (OPTIONAL) Tabix index for TRUTH_VCF. Made if not provided.
- *EVALUATION_REGIONS_BED*: BED of regions to evaluate in
- *STRATIFICATION_ARCHIVE*: (OPTIONAL) tar.gz of a GIAB-style stratification folder (root TSV plus its referenced BED files) to break the Aardvark results down by
- *RESTRICT_REGIONS_BED*: (OPTIONAL) Additional BED to restrict comparison against TRUTH_VCF to
- *OUTPUT_GAF*: Should a GAF file with the aligned reads be saved for each run? When both runs map the same way there is only one set of alignments, so both outputs are the same file. Default is 'false'.
- *OUTPUT_BAM*: Should the merged BAM be saved for each run? Default is 'false'.
- *PAIRED_READS*: Are the reads paired? Default is 'true'.
- *INTERLEAVED_READS*: Are paired reads interleaved in a single FASTQ? Only meaningful when PAIRED_READS is true and there is a single input FASTQ. Default is 'false'.
- *READS_PER_CHUNK*: Number of reads contained in each mapping chunk. Default 20 million.
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
- *DV_MODEL_TYPE*: Type of DeepVariant model to use. Can be WGS (default), WES, PACBIO, ONT_R104, or HYBRID_PACBIO_ILLUMINA.
- *DV_MODEL_META*: .meta file for a custom DeepVariant calling model
- *DV_MODEL_INDEX*: .index file for a custom DeepVariant calling model
- *DV_MODEL_DATA*: .data-00000-of-00001 file for a custom DeepVariant calling model
- *DV_MODEL_FILES*: Array of all files in the root directory of the DV model, if not using DV_MODEL_META/DV_MODEL_INDEX/DV_MODEL_DATA format
- *DV_MODEL_VARIABLES_FILES*: Array of files that need to go in a 'variables' subdirectory for a DV model
- *DV_KEEP_LEGACY_AC*: Should DV use the legacy allele counter behavior? If unspecified this is not done, unless set in the model. Might want to be on for short reads.
- *DV_NORM_READS*: Should DV normalize reads itself? If unspecified this is not done, unless set in the model.
- *OTHER_MAKEEXAMPLES_ARG*: Additional arguments for the make_examples step of DeepVariant
- *DV_USE_GPUS*: Should DeepVariant use GPUs for calling variants? Default is 'true'.
- *DV_NO_GPU_DOCKER*: Container image to use when running DeepVariant for steps that don't benefit from GPUs. Must be DeepVariant 1.8+.
- *DV_GPU_DOCKER*: Container image to use when running DeepVariant for steps that benefit from GPUs. Must be DeepVariant 1.8+.
- *SPLIT_READ_CORES*: Number of cores to use when splitting the reads into chunks. Default is 8.
- *SPLIT_READ_MEM*: Memory, in GB, to use when splitting the reads into chunks. Default is 50.
- *MAP_CORES*: Number of cores to use when mapping the reads. Default is 16.
- *MAP_MEM*: Memory, in GB, to use when mapping the reads. Default is 120.
- *HAPLOTYPE_SAMPLING*: Whether or not to use haplotype sampling before running giraffe. The sampled graph and its indexes count as indexes, so they are made once unless CANDIDATE_SEPARATE_INDEXES is set. Default is 'true'.
- *SET_REFERENCE*: (OPTIONAL) Name of the single reference to keep for haplotype sampling.
- *INDEX_MINIMIZER_WEIGHTED*: Whether to use weighted minimizer indexing. (Default: true)
- *INDEX_MINIMIZER_MEM*: Memory, in GB, to use when making the minimizer index. (Default: 320 if weighted, 120 otherwise)
- *KMER_COUNTING_MEM*: Memory, in GB, to use when counting kmers. (Default: 64)
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling, and giraffe distance index). (Default: 120)
- *BAM_PREPROCESS_MEM*: Memory, in GB, to use when preprocessing BAMs (left-shifting and preparing realignment targets). Default is 20.
- *REALIGN_MEM*: Memory, in GB, to use for Abra indel realignment. Default is 40 or MAP_MEM, whichever is lower.
- *CALL_CORES*: Number of cores to use when calling variants. Default is 8.
- *CALL_MEM*: Memory, in GB, to use when calling variants. Default is 50.
- *MAKE_EXAMPLES_CORES*: Number of cores to use when making DeepVariant examples. Default is CALL_CORES.
- *MAKE_EXAMPLES_MEM*: Memory, in GB, to use when making DeepVariant examples. Default is CALL_MEM.
- *EVAL_CORES*: Number of cores to use when evaluating variant calls. Default is 16.
- *EVAL_MEM*: Memory, in GB, to use when evaluating variant calls. Default is 30.

Related
topics: [read realignment](#Read-realignment), [reference prefix removal](#Reference-prefix-removal), [CRAM input](#CRAM-input), [reads chunking](#Reads-chunking), [path list](#Path-list), [single-end reads](#Single-end-reads), [interleaved reads](#Interleaved-reads), [HPRC pangenomes](#HPRC-pangenomes).

[Test locally](#testing-locally) with:

```sh
miniwdl run --as-me workflows/acceptance_test.wdl -i params/acceptance_test.json
```

### GAF to sorted GAM workflow

Currently, only GAM file can be sorted and indexed, for example to extract and subgraph and visualize, or use with
the [sequenceTubeMap](https://github.com/vgteam/sequenceTubeMap).
This workflow converts reads aligned to a pangenome in a GAF file to a sorted and indexed GAM file.

- workflow file: [workflows/sort_graph_aligned_reads.wdl](workflows/sort_graph_aligned_reads.wdl)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/sortGraphAlignedReads:master?tab=info)

Parameters:

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

Parameters:

- *GBZ_FILE*: Path to .gbz index file
- *INPUT_READ_FILE_FIRST*: Input sample 1st read pair fastq.gz or fastq
- *INPUT_READ_FILE_SECOND*: Input sample 2nd read pair fastq.gz or fastq
- *HAPL_FILE*: Path to .hapl file
- *DIST_FILE*: Path to .dist file
- *R_INDEX_FILE*: Path to .ri file
- *KFF_FILE*: Path to .kff file
- *OUTPUT_NAME_PREFIX*: Name of the output file (Default: haplotype_sampled_graph)
- *KMER_LENGTH*: Size of kmer using for sampling (Up to 31) (Default: 29)
- *CORES*: Number of cores to use with commands. (Default: 16)
- *KMER_COUNTING_MEM*: Memory, in GB, to use when counting kmers. (Default: 64)
- *HAPLOTYPE_INDEXING_MEM*: Memory, in GB, to use for haplotype sampling indexing tasks (distance index, r-index, haplotype index, sampling). (Default: 120)
- *WINDOW_LENGTH*: Window length used for building the minimizer index for sampling haplotypes. (Default: 11)
- *SUBCHAIN_LENGTH*: Target length (in bp) for subchains. (Default: 10000)
- *HAPLOTYPE_NUMBER*: Number of generated synthetic haplotypes. (Default: 4)
- *PRESENT_DISCOUNT*: Multiplicative factor for discounting scores for present kmers. (Default: 0.9)
- *HET_ADJUST*: Additive term for adjusting scores for heterozygous kmers. (Default: 0.05)
- *ABSENT_SCORE*: Score for absent kmers. (Default: 0.8)
- *INCLUDE_REFERENCE*: Include reference paths and generic paths from the full graph in the sampled graph. (Default: true)
- *SET_REFERENCE*: Name of single reference to include in sampled graph. (Default: all references)
- *DIPLOID*: Activate diploid sampling. (Default: true)
- *OUTPUT_HAPL*: Whether or not to output the .hapl haplotype information file computed before sampling, so a caller mapping and calling against the same pangenome can reuse it instead of recomputing it. Default is 'false'.
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

### Internal subworkflows

[workflows/internal/](workflows/internal) holds subworkflows that exist only to be called by the workflows above. They
are steps that several workflows need to do the same way, pulled out so there is one copy of each, and they are not
meant to be run on their own. Their parameters are documented in their own `parameter_meta` sections rather than here,
and they can change without notice.

- [workflows/internal/split_reads.wdl](workflows/internal/split_reads.wdl) (`SplitReads`): get a sample's reads as FASTQ,
  converting a CRAM or a BAM if that is what arrived, and split them into chunks to map in parallel.
- [workflows/internal/prepare_reference.wdl](workflows/internal/prepare_reference.wdl) (`PrepareReference`): work out
  which contigs to work on and get a FASTA reference for them, with its `.fai` and `.dict`, extracting them from the
  graph if they weren't provided.
- [workflows/internal/surject.wdl](workflows/internal/surject.wdl) (`Surject`): project GAF chunks onto the reference
  paths with `vg surject` and merge them into one sorted BAM.
- [workflows/internal/index_for_giraffe.wdl](workflows/internal/index_for_giraffe.wdl) (`IndexForGiraffe`): produce the
  GBZ, distance, minimizer and zipcodes indexes `vg giraffe` needs, building whatever wasn't provided and haplotype
  sampling the graph if asked.

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

Continuous integration runs `miniwdl check` on every workflow, and separately runs:

```sh
python3 scripts/lint_wdl_docs.py
```

which checks that every workflow input is described in the workflow's `parameter_meta` section, and that every workflow
meant to be run has a section in this README listing all of its parameters.
[Internal subworkflows](#internal-subworkflows) still need `parameter_meta`, but are not expected to have a README
section. Files that have never satisfied any of this are listed
in [scripts/doc_lint_exemptions.txt](scripts/doc_lint_exemptions.txt); new workflows are expected not to need an entry
there.

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
