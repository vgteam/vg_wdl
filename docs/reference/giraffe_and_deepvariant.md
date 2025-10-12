## Reference: Inputs, outputs, and defaults for `GiraffeDeepVariant`

### Required inputs
- `GBZ_FILE` (File): GBZ graph.
- `DIST_FILE` (File): Distance index.
- `MIN_FILE` (File): Minimizer index.
- `SAMPLE_NAME` (String): Sample identifier.
- Reads: one of
  - FASTQs: `INPUT_READ_FILE_1` (File), `INPUT_READ_FILE_2` (File, optional if single‑end; set `PAIRED_READS=false`), or
  - CRAM: `INPUT_CRAM_FILE` (File), plus `CRAM_REF` (File) and `CRAM_REF_INDEX` (File)

### Optional inputs (type; default)
- `ZIPCODES_FILE` (File?; none)
- `OUTPUT_GAF` (Boolean; true)
- `OUTPUT_SINGLE_BAM` (Boolean; false)
- `OUTPUT_CALLING_BAMS` (Boolean; `!OUTPUT_SINGLE_BAM`)
- `OUTPUT_UNMAPPED_BAM` (Boolean; false)
- `PAIRED_READS` (Boolean; true)
- `READS_PER_CHUNK` (Int; 20000000)
- `CONTIGS` (Array[String]+?; none)
- `PATH_LIST_FILE` (File?; none)
- `REFERENCE_PREFIX` (String; "")
- `REFERENCE_FILE` (File?; none)
- `REFERENCE_INDEX_FILE` (File?; none)
- `REFERENCE_DICT_FILE` (File?; none)
- `HAPLOID_CONTIGS` (Array[String]?; none)
- `PAR_REGIONS_BED_FILE` (File?; none)
- `PRUNE_LOW_COMPLEXITY` (Boolean; true)
- `LEFTALIGN_BAM` (Boolean; true)
- `REALIGN_INDELS` (Boolean; true)
- `REALIGNMENT_EXPANSION_BASES` (Int; 160)
- `MIN_MAPQ` (Int?; model default)
- `MAX_FRAGMENT_LENGTH` (Int; 3000)
- `GIRAFFE_PRESET` (String; "default")
- `GIRAFFE_OPTIONS` (String; "")
- `TRUTH_VCF` (File?; none)
- `TRUTH_VCF_INDEX` (File?; none)
- `EVALUATION_REGIONS_BED` (File?; none)
- `RESTRICT_REGIONS_BED` (File?; none)
- `TARGET_REGION` (String?; none)
- `RUN_STANDALONE_VCFEVAL` (Boolean; true)
- `DV_MODEL_TYPE` (String; "WGS")
- `DV_MODEL_META` (File?; none)
- `DV_MODEL_INDEX` (File?; none)
- `DV_MODEL_DATA` (File?; none)
- `DV_MODEL_FILES` (Array[File]?; none)
- `DV_MODEL_VARIABLES_FILES` (Array[File]?; none)
- `DV_KEEP_LEGACY_AC` (Boolean?; model default)
- `DV_NORM_READS` (Boolean?; model default)
- `OTHER_MAKEEXAMPLES_ARG` (String; "")
- `DV_NO_GPU_DOCKER` (String?; none)
- `DV_GPU_DOCKER` (String?; none)
- `SPLIT_READ_CORES` (Int; 8)
- `MAP_CORES` (Int; 16)
- `MAP_MEM` (Int; 120)  // GB
- `CALL_CORES` (Int; 8)
- `CALL_MEM` (Int; 50)   // GB
- `MAKE_EXAMPLES_CORES` (Int; `CALL_CORES`)
- `MAKE_EXAMPLES_MEM` (Int; `CALL_MEM`)
- `EVAL_MEM` (Int; 60)   // GB
- `VG_DOCKER` (String; "quay.io/vgteam/vg:v1.64.0")
- `VG_GIRAFFE_DOCKER` (String?; none)
- `VG_SURJECT_DOCKER` (String?; none)

### Outputs
- `output_vcf` (File): VCF with variant calls
- `output_vcf_index` (File): Tabix index for VCF
- `output_gvcf` (File): gVCF
- `output_gvcf_index` (File): Tabix index for gVCF
- `output_gaf` (File?; present if `OUTPUT_GAF=true`): GAF alignments
- `output_bam` (File?): Merged BAM (if enabled)
- `output_bam_index` (File?): BAM index (if enabled)
- `output_calling_bams` (Array[File]?): Per‑contig calling BAMs (if enabled)
- `output_calling_bam_indexes` (Array[File]?): Indexes for calling BAMs
- `output_unmapped_bam` (File?): Unmapped reads BAM (if enabled)
- `output_vcfeval_evaluation_archive` (File?): vcfeval results archive (if truth provided)
- `output_happy_evaluation_archive` (File?): hap.py results archive (if truth provided)

### Example input JSON (FASTQs)

```json
{
  "GiraffeDeepVariant.INPUT_READ_FILE_1": "tests/small_sim_graph/reads_1.fastq.gz",
  "GiraffeDeepVariant.INPUT_READ_FILE_2": "tests/small_sim_graph/reads_2.fastq.gz",
  "GiraffeDeepVariant.GBZ_FILE": "tests/small_sim_graph/graph.gbz",
  "GiraffeDeepVariant.SAMPLE_NAME": "s0",
  "GiraffeDeepVariant.MIN_FILE": "tests/small_sim_graph/graph.min",
  "GiraffeDeepVariant.DIST_FILE": "tests/small_sim_graph/graph.dist",
  "GiraffeDeepVariant.DV_NO_GPU_DOCKER": "google/deepvariant:1.8.0",
  "GiraffeDeepVariant.VG_DOCKER": "quay.io/vgteam/vg:v1.64.0"
}
```


