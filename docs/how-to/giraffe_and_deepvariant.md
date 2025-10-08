## How‑Tos: Giraffe + DeepVariant

Targeted guides for common tasks and configurations.

### Use FASTQs (paired‑end) vs. CRAM

- **FASTQs (paired)**: Provide `INPUT_READ_FILE_1`, `INPUT_READ_FILE_2`, and set `PAIRED_READS=true` (default).
- **FASTQ (single‑end)**: Provide only `INPUT_READ_FILE_1` and set `PAIRED_READS=false`.
- **CRAM**: Provide `INPUT_CRAM_FILE` plus `CRAM_REF` and `CRAM_REF_INDEX` (reference FASTA and `.fai`). Omit FASTQ inputs.

### Choose graph contigs / reference naming

- To restrict calling to particular paths, either:
  - Set `CONTIGS` to an array of path names, or
  - Provide `PATH_LIST_FILE` with one path name per line.
- If not set, the workflow extracts major path names from the `GBZ_FILE` and filters out decoys.
- Use `REFERENCE_PREFIX` to strip a prefix from graph path names when creating the reference and surjecting BAM (e.g., to match `chr` prefixes).

### Toggle outputs

- `OUTPUT_GAF` (default `true`): write GAF alignments.
- `OUTPUT_SINGLE_BAM` (default `false`): write a single merged BAM. When `true`, `OUTPUT_CALLING_BAMS` defaults to `false`.
- `OUTPUT_CALLING_BAMS`: write per‑contig BAMs used for calling.
- `OUTPUT_UNMAPPED_BAM` (default `false`).

### CPU vs GPU DeepVariant

- CPU image: set `DV_NO_GPU_DOCKER`, for example `google/deepvariant:1.8.0`.
- GPU image: set `DV_GPU_DOCKER`, for example `google/deepvariant:1.8.0-gpu`.
- Don’t set both simultaneously; provide whichever matches your environment.
- Select model in `DV_MODEL_TYPE`: one of `WGS` (default), `WES`, `PACBIO`, `ONT_R104`, `HYBRID_PACBIO_ILLUMINA`.
- To use a custom model, supply either `DV_MODEL_META`/`DV_MODEL_INDEX`/`DV_MODEL_DATA` or `DV_MODEL_FILES` (+ optional `DV_MODEL_VARIABLES_FILES`).

### Resource tuning

- Giraffe mapping: `MAP_CORES` (default 16), `MAP_MEM` GB (default 120), `READS_PER_CHUNK` (default 20,000,000).
- DeepVariant: `CALL_CORES` (default 8), `CALL_MEM` GB (default 50), `MAKE_EXAMPLES_*` (default to CALL_*), `EVAL_MEM` GB (default 60).
- Indel realignment limits memory internally with `REALIGN_MEM = min(MAP_MEM, 40)`.

### CRAM specifics

- When `INPUT_CRAM_FILE` is set, you must provide `CRAM_REF` and `CRAM_REF_INDEX` that correspond to the CRAM’s reference.
- If your graph does not include all reference bases, provide `REFERENCE_FILE`/`REFERENCE_INDEX_FILE`/`REFERENCE_DICT_FILE` explicitly.

### Add truth VCF evaluation

- Provide `TRUTH_VCF` and `TRUTH_VCF_INDEX` (tabix) and optionally region controls:
  - `EVALUATION_REGIONS_BED` and/or `RESTRICT_REGIONS_BED`, `TARGET_REGION`.
- `RUN_STANDALONE_VCFEVAL` runs vcfeval separately from hap.py; keep `true` unless you observe incompatibilities.

### Run with miniwdl (local)

```bash
miniwdl run workflows/giraffe_and_deepvariant.wdl \
  -i params/giraffe_and_deepvariant.json
```


