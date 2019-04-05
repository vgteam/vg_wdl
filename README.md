vg\_wdl
---------------
Eric T Dawson, Mike Lin and Charles Markello
MIT License, 2018

Workflow Description Language (WDL) scripts for common vg workflows

## Overview

## Docker Containers
WDL needs the runtime Docker image to be present on Dockerhub.  
VG images are available in [quay](https://quay.io/repository/vgteam/vg?tab=tags)
and selected images are available in [ the variantgraphs Dockerhub](https://cloud.docker.com/u/variantgraphs/repository/docker/variantgraphs/vg),  
and can be pulled with:  
```
    docker pull variantgraphs/vg  
```

Specific tags can be specified like so:  
```
# Get the vg 1.3.1 release  
docker pull variantgraphs/vg:1.3.1
```

## Usage

### NIH Biowulf Usage
#### Tiny chromosome 21 example
cd into a directory with at least 200 GB of allocated Disk space
```
cd /data/$USER
```

Launch an interactive session on Biowulf:
```
sinteractive -c6 --mem=100g --gres=lscratch:50 --time=20:00:00
```

Load requisite Biowulf modules:
```
module load cromwell/36 git
```

Clone the github repo and create a work directory for running the wdl workflow:
```
git clone https://github.com/vgteam/vg_wdl.git
mkdir -p vg_wdl_testrun && cd vg_wdl_testrun
```

Copy requisite wdl script files into work directory:
```
cp ../vg_wdl/workflows/vg_multi_map_call.wdl ../vg_wdl/workflows/vg_multi_map_call.inputs_tiny_chr21.biowulf_example.json ../vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf .
```

Create an input directory and download example input files:
```
mkdir -p workflow_inputs
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG002_chr21_1.tiny.fastq.gz -O workflow_inputs/HG002_chr21_1.tiny.fastq.gz
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG002_chr21_2.tiny.fastq.gz -O workflow_inputs/HG002_chr21_2.tiny.fastq.gz
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/path_list_21.txt -O workflow_inputs/path_list_21.txt
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/path_list_length_21.txt -O workflow_inputs/path_list_length_21.txt
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_chr21_t289.xg -O workflow_inputs/snp1kg_maf0.01_chr21_t289.xg
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_chr21_t289.gcsa -O workflow_inputs/snp1kg_maf0.01_chr21_t289.gcsa
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_chr21_t289.gcsa.lcp -O workflow_inputs/snp1kg_maf0.01_chr21_t289.gcsa.lcp
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_chr21_t289.gbwt -O workflow_inputs/snp1kg_maf0.01_chr21_t289.gbwt
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa -O workflow_inputs/hs37d5.fa
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa.fai -O workflow_inputs/hs37d5.fa.fai
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.dict -O workflow_inputs/hs37d5.dict
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/snpEff_v4_3_GRCh37.75.zip -O workflow_inputs/snpEff_v4_3_GRCh37.75.zip
```

Run the workflow:
```
java -Dconfig.file=custom_biowulf_cromwell_singularity.conf -jar ${CROMWELL_JAR} run -i vg_multi_map_call.inputs_tiny_chr21.biowulf_example.json vg_multi_map_call.wdl
```

#### Whole genome example
cd into a directory with at least 200 GB of allocated Disk space
```
cd /data/$USER
```

Launch an interactive session on Biowulf:
```
sinteractive -c6 --mem=100g --gres=lscratch:50 --time=20:00:00
```

Load requisite Biowulf modules:
```
module load cromwell/36 git
```

Clone the github repo and create a work directory for running the wdl workflow:
```
git clone https://github.com/vgteam/vg_wdl.git
mkdir -p vg_wdl_testrun && cd vg_wdl_testrun
```

Copy requisite wdl script files into work directory:
```
cp ../vg_wdl/workflows/vg_multi_map_call.wdl ../vg_wdl/workflows/vg_multi_map_call.inputs_whole_genome.biowulf_example.json ../vg_wdl/workflows/custom_biowulf_cromwell_singularity.conf .
```
Create an input directory and download example input files:
```
mkdir -p workflow_inputs
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG002_1.sorted.fastq.gz -O workflow_inputs/HG002_1.sorted.fastq.gz
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG002_2.sorted.fastq.gz -O workflow_inputs/HG002_2.sorted.fastq.gz
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/path_list_whole_genome.txt -O workflow_inputs/path_list_whole_genome.txt
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/path_list_length_whole_genome.txt -O workflow_inputs/path_list_length_whole_genome.txt
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/snp1kg_maf0.01_decoys_minaf_0.01.xg -O workflow_inputs/snp1kg_maf0.01_decoys_minaf_0.01.xg
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/snp1kg_maf0.01_decoys_minaf_0.01.gcsa -O workflow_inputs/snp1kg_maf0.01_decoys_minaf_0.01.gcsa
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/snp1kg_maf0.01_decoys_minaf_0.01.gcsa.lcp -O workflow_inputs/snp1kg_maf0.01_decoys_minaf_0.01.gcsa.lcp
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/snp1kg_maf0.01_decoys_minaf_0.01.gbwt -O workflow_inputs/snp1kg_maf0.01_decoys_minaf_0.01.gbwt
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa -O workflow_inputs/hs37d5.fa
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa.fai -O workflow_inputs/hs37d5.fa.fai
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.dict -O workflow_inputs/hs37d5.dict
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/snpEff_v4_3_GRCh37.75.zip -O workflow_inputs/snpEff_v4_3_GRCh37.75.zip
```

Run the workflow:
```
java -Dconfig.file=custom_biowulf_cromwell_singularity.conf -jar ${CROMWELL_JAR} run -i vg_multi_map_call.inputs_whole_genome.biowulf_example.json vg_multi_map_call.wdl
```

When running the workflow on new data, be sure to set the values of `vgPipeline.INPUT_READ_FILE_1` and `vgPipeline.INPUT_READ_FILE_2` in the input JSON file (ex: `vg_multi_map_call.inputs_tiny_chr21.biowulf_example.json` to the full paths of your respective read pair FASTQ files.

NOTE IF USING DRAGEN MODULE:
* Set `vgPipeline.UDPBINFO_PATH` in the input JSON file to a path in a Udp directory that also includes your Biowulf username. For example set `vgPipeline.UDPBINFO_PATH` to `Udpwork/usr/markellocj` if your biowulf username is `markellocj`.
* Set `vgPipeline.HELIX_USERNAME` in the input JSON file to your Biowulf username. For example set `vgPipeline.HELIX_USERNAME` to `markellocj` if your biowulf username is `markellocj`.

#### Pipeline options
Mapping algorithm options:
* To run the workflow using VG MAP mapping algorithm, in the input JSON file set `vgPipeline.RUN_VGMPMAP_ALGORITHM` to `false`.
* To run the workflow using VG MPMAP mapping algorithm, in the input JSON file set `vgPipeline.RUN_VGMPMAP_ALGORITHM` to `true`.

Caller algorithm options:
* To run the workflow using the VG CALL calling algorithm, in the input JSON file set `vgPipeline.RUN_LINEAR_CALLER` to `false`.
* To run the workflow using GATKs HaplotypeCaller calling algorithm, in the input JSON file set `vgPipeline.RUN_LINEAR_CALLER` to `true` AND set `vgPipeline.RUN_DRAGEN_CALLER` to `false`.
* To run the workflow using the DRAGEN FPGA module for calling genotypes, in the input JSON file set `vgPipeline.RUN_LINEAR_CALLER` to `true` AND set `vgPipeline.RUN_DRAGEN_CALLER` to `true`.

## Examples

## Availability

## Contributing, Help, Bugs and Requests
Please open an Issue on [github](https://github.com/vgteam/vg_wdl) for help, bug reports, or feature requests.
When doing so, please remember that vg\_wdl is open-source software made by a community of developers. Please
be considerate and support a positive environment.
