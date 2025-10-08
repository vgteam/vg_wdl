## Giraffe + DeepVariant: Hands‑on Tutorial

This tutorial walks you through running the `GiraffeDeepVariant` workflow end‑to‑end on the included small test dataset. You’ll map reads to a pangenome with vg Giraffe and call small variants with DeepVariant.

### Prerequisites
- Docker or Singularity available on your system
- miniwdl installed

### Test data and key files
- Inputs JSON: `params/giraffe_and_deepvariant.json`
- Workflow: `workflows/giraffe_and_deepvariant.wdl`
- Small test dataset: `tests/small_sim_graph/` (graph indices and reads)

The example inputs reference:
- `GBZ_FILE`: `tests/small_sim_graph/graph.gbz`
- `MIN_FILE`: `tests/small_sim_graph/graph.min`
- `DIST_FILE`: `tests/small_sim_graph/graph.dist`
- `INPUT_READ_FILE_1`: `tests/small_sim_graph/reads_1.fastq.gz`
- `INPUT_READ_FILE_2`: `tests/small_sim_graph/reads_2.fastq.gz`

### Run with miniwdl

```bash
cd /Users/seeskand/PycharmProjects/vg_wdl_fork
miniwdl run workflows/giraffe_and_deepvariant.wdl \
  -i params/giraffe_and_deepvariant.json
```

Notes:
- The default params file sets `GiraffeDeepVariant.DV_GPU_DOCKER` to a GPU-enabled DeepVariant image. If you don’t have a GPU, set `GiraffeDeepVariant.DV_NO_GPU_DOCKER` to `google/deepvariant:1.8.0` and unset `DV_GPU_DOCKER` in your input JSON.

Tips:
- Add `--writeLogs $(pwd)/.toil-logs` and/or `--logDebug` for detailed logs.
- To resume a failed run: `toil restart file:$(pwd)/gdv-jobstore`.

### What you get (outputs)
The workflow produces, at a minimum:
- `output_vcf` and `output_vcf_index`: small-variant VCF + index
- `output_gvcf` and `output_gvcf_index`: gVCF + index
- `output_bam` and `output_bam_index` (if enabled via output toggles)
- `output_gaf` from Giraffe mapping if `OUTPUT_GAF=true`
- Optional evaluation archives if truth data are provided

miniwdl writes a run directory containing `outputs.json` with paths to all products.

### Next steps
- See the How‑Tos to customize inputs (FASTQs vs CRAM, output toggles, tuning CPU/memory, GPU vs CPU DeepVariant).
- See the Reference for the full list of inputs/outputs with types and defaults.


