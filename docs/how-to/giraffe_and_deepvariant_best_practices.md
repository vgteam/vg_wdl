## Best Practices: Giraffe + DeepVariant

Use this simple, production‑oriented recipe to run `GiraffeDeepVariant` on your own data. Replace only the input read paths; keep the rest as recommended defaults unless you know you need to change them.

### 1) Download required graph/indexes

These files are large; use a server or a reliable network. Download into `graphs/`:

```bash
mkdir -p graphs && cd graphs

# GBZ (graph)
curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz

# Minimizer index
curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.min

# Distance index
curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.dist

```

Files used above:
- [HPRC v1.1 GBZ](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz)
- [HPRC v1.1 MIN](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.min)
- [HPRC v1.1 DIST](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.dist)

### 2) Prepare inputs JSON (template)

Save the following as `my_inputs.json`, then edit only the FASTQ/CRAM inputs and `SAMPLE_NAME`:

```json
{
  "GiraffeDeepVariant.INPUT_READ_FILE_1": "/path/to/reads_R1.fastq.gz",
  "GiraffeDeepVariant.INPUT_READ_FILE_2": "/path/to/reads_R2.fastq.gz",

  "GiraffeDeepVariant.GBZ_FILE": "graphs/hprc-v1.1-mc-grch38.gbz",
  "GiraffeDeepVariant.MIN_FILE": "graphs/hprc-v1.1-mc-grch38.min",
  "GiraffeDeepVariant.DIST_FILE": "graphs/hprc-v1.1-mc-grch38.dist",

  "GiraffeDeepVariant.REFERENCE_PREFIX": "GRCh38#0#",

  "GiraffeDeepVariant.SAMPLE_NAME": "sample_name",

  "GiraffeDeepVariant.OUTPUT_GAF": true,
  "GiraffeDeepVariant.OUTPUT_SINGLE_BAM": false,

  "GiraffeDeepVariant.DV_MODEL_TYPE": "WGS",
  "GiraffeDeepVariant.DV_GPU_DOCKER": "google/deepvariant:1.8.0-gpu",
  "GiraffeDeepVariant.VG_DOCKER": "quay.io/vgteam/vg:v1.68.0",

  "GiraffeDeepVariant.MAP_CORES": 16,
  "GiraffeDeepVariant.MAP_MEM": 120,
  "GiraffeDeepVariant.CALL_CORES": 8,
  "GiraffeDeepVariant.CALL_MEM": 50
}
```

Notes:
- For single‑end reads, omit `INPUT_READ_FILE_2` and add `"GiraffeDeepVariant.PAIRED_READS": false`.
- For CRAM input, set `GiraffeDeepVariant.INPUT_CRAM_FILE`, `CRAM_REF`, and `CRAM_REF_INDEX`, and remove FASTQ inputs. 
- Set `"GiraffeDeepVariant.REFERENCE_PREFIX"` to `"GRCh38#0#"` for GRCh38-based graphs, and to `"CHM13#0#"` for CHM13-based graphs.

### Worked example: HG002 public reads (40x PCR‑free)

If you want to run the pipeline for a public sample, download the public FASTQs into `reads/`:

```bash
mkdir -p reads && cd reads
curl -O https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG002.novaseq.pcr-free.40x.R1.fastq.gz
curl -O https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG002.novaseq.pcr-free.40x.R2.fastq.gz
cd -
```



Create `hg002_inputs.json`:

```json
{
  "GiraffeDeepVariant.INPUT_READ_FILE_1": "reads/HG002.novaseq.pcr-free.40x.R1.fastq.gz",
  "GiraffeDeepVariant.INPUT_READ_FILE_2": "reads/HG002.novaseq.pcr-free.40x.R2.fastq.gz",

  "GiraffeDeepVariant.GBZ_FILE": "graphs/hprc-v1.1-mc-grch38.gbz",
  "GiraffeDeepVariant.MIN_FILE": "graphs/hprc-v1.1-mc-grch38.min",
  "GiraffeDeepVariant.DIST_FILE": "graphs/hprc-v1.1-mc-grch38.dist",

  "GiraffeDeepVariant.REFERENCE_PREFIX": "GRCh38#0#",

  "GiraffeDeepVariant.SAMPLE_NAME": "HG002",

  "GiraffeDeepVariant.OUTPUT_GAF": true,
  "GiraffeDeepVariant.OUTPUT_SINGLE_BAM": false,

  "GiraffeDeepVariant.DV_MODEL_TYPE": "WGS",
  "GiraffeDeepVariant.DV_GPU_DOCKER": "google/deepvariant:1.8.0-gpu",
  "GiraffeDeepVariant.VG_DOCKER": "quay.io/vgteam/vg:v1.68.0",

  "GiraffeDeepVariant.MAP_CORES": 16,
  "GiraffeDeepVariant.MAP_MEM": 120,
  "GiraffeDeepVariant.CALL_CORES": 8,
  "GiraffeDeepVariant.CALL_MEM": 50
}
```

Run:

```bash
miniwdl run workflows/giraffe_and_deepvariant.wdl \
  -i hg002_inputs.json
```

### 3) Run (miniwdl)

```bash
miniwdl run workflows/giraffe_and_deepvariant.wdl \
  -i my_inputs.json
```

Outputs are listed in the run directory’s `outputs.json`.

### Recommended resources

- GPU strongly recommended for DeepVariant: use a GPU‑enabled container (`DV_GPU_DOCKER`) when available.
- Memory/CPU guidelines:
  - Mapping: `MAP_CORES=16`, `MAP_MEM=120` (GB) are good defaults for WGS scale.
  - Calling: `CALL_CORES=8`, `CALL_MEM=50` (GB) are typical starting points.

### Minimal changes for your data

- Change only:
  - `INPUT_READ_FILE_1` (and `INPUT_READ_FILE_2` if paired)
  - `SAMPLE_NAME`
- Keep graph/index paths pointing to your chosen reference graph.
- The default params file sets `GiraffeDeepVariant.DV_GPU_DOCKER` to a GPU-enabled DeepVariant image. If you don’t have a GPU, set `GiraffeDeepVariant.DV_NO_GPU_DOCKER` to `google/deepvariant:1.8.0` and unset `DV_GPU_DOCKER` in your input JSON.


### Notes for running on a Slurm system using toil-wdl-runner
Use the folowing command: 
```bash
toil-wdl-runner --jobStore ./big_store --batchSystem slurm --caching false --batchLogsDir ./logs workflows/giraffe_and_deepvariant.wdl hg002_inputs.json -o slurm_run -m slurm_run.
json
```
- If you are using gpu activated node, you need to use export TOIL_SLURM_ARGS="--time=12:00:00 --partition=gpu" (only if the server only give gpu access that way)

