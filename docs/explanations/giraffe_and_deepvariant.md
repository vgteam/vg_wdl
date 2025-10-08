## Explanations: How the Giraffe + DeepVariant workflow works

### Overview
The workflow maps sequencing reads to a pangenome graph using vg Giraffe and then calls small variants with DeepVariant. It builds a linear reference FASTA from the graph paths you specify to ensure coordinate consistency between mapping, surjection to BAM, and variant calling.

### Graph inputs and path selection
- `GBZ_FILE` contains the GBZ graph with path names (e.g., chromosomes). `DIST_FILE` and `MIN_FILE` are additional indexes required by Giraffe.
- If you do not provide `CONTIGS` or `PATH_LIST_FILE`, the workflow extracts a list of “major” paths from the GBZ and filters out decoys/unlocalized contigs to keep the scatter practical and evaluation meaningful.
- `REFERENCE_PREFIX` lets you strip a prefix from path names (e.g., remove `chr` or add compatibility) when generating the FASTA and when surjecting BAM so that names match downstream tooling.

### Reference generation and indexing
Unless `REFERENCE_FILE` is provided, the workflow extracts a FASTA reference from the selected graph paths. It then indexes it (`.fai`) and builds a sequence dictionary (`.dict`). Providing your own reference can be necessary if the graph lacks bases for certain paths.

### Mapping and BAM preparation
- Giraffe aligns reads to the graph using presets (`GIRAFFE_PRESET`: `default`, `fast`, `hifi`, or `r10`) and optional extra flags (`GIRAFFE_OPTIONS`).
- Reads are then surjected to the linear reference and optionally post‑processed:
  - Low‑complexity pruning during surjection (`PRUNE_LOW_COMPLEXITY`)
  - Left‑align indels (`LEFTALIGN_BAM`)
  - Targeted indel realignment (`REALIGN_INDELS`, with window size `REALIGNMENT_EXPANSION_BASES`)
- Output toggles control whether a single merged BAM, per‑contig calling BAMs, or an unmapped BAM are emitted.

### DeepVariant calling
- Model: `DV_MODEL_TYPE` selects among tuned models; custom models can be supplied via TensorFlow checkpoint files (`DV_MODEL_META`, `DV_MODEL_INDEX`, `DV_MODEL_DATA`, or `DV_MODEL_FILES`). Optional `DV_MODEL_VARIABLES_FILES` populate a `variables/` directory expected by some models.
- GPU vs CPU: set one of `DV_GPU_DOCKER` (GPU‑accelerated steps) or `DV_NO_GPU_DOCKER` (CPU). Only steps that benefit from GPUs will use the GPU image.
- MAPQ filtering (`MIN_MAPQ`) can be set or left to DeepVariant defaults for the model.

### Evaluation (optional)
If ground truth is supplied (`TRUTH_VCF`, `TRUTH_VCF_INDEX`), the workflow can evaluate calls using hap.py, optionally restricting regions (`EVALUATION_REGIONS_BED`, `RESTRICT_REGIONS_BED`, `TARGET_REGION`). `RUN_STANDALONE_VCFEVAL` also runs vcfeval separately, though some DeepVariant VCFs may not be compatible in older versions.

### Special samples/regions
- `HAPLOID_CONTIGS` and `PAR_REGIONS_BED_FILE` control sex‑chromosome ploidy handling; note incompatibility with DeepVariant 1.5 for these inputs.

### Resources and scalability
- Mapping is memory‑intensive (`MAP_MEM`, default 120 GB) and can use many cores (`MAP_CORES`, default 16). DeepVariant has separate CPU/memory controls for make‑examples and call‑variants steps.
- The workflow scatters across contigs/paths; choosing fewer, larger paths balances parallelism with overhead.


