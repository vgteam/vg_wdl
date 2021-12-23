vg\_wdl
---------------
Eric T Dawson, Mike Lin and Charles Markello, Jean Monlong, Adam Novak
MIT License, 2022

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

Install miniwdl in a python 3 virtual environment
```
git clone https://github.com/chanzuckerberg/miniwdl.git
virtualenv miniwdl_venv
source miniwdl_venv/bin/activate
pip3 install ./miniwdl
deactivate
```
Download `.wdl` workflow file and `.json` input file
```
wget https://github.com/vgteam/vg_wdl/raw/master/workflows/vg_multi_map_call.wdl
wget https://github.com/vgteam/vg_wdl/raw/master/params/vg_multi_map_call.inputs_tiny.http_url.json
```
Activate the miniwdl virtual environment and run the example workflow
```
source miniwdl_venv/bin/activate
miniwdl cromwell vg_multi_map_call.wdl -i vg_multi_map_call.inputs_tiny.http_url.json
```
To modify the input parameters, edit the input `.json` with the necessary changes.

## Examples

## Availability

Workflow for processing single sample datasets:
- [workflow file](https://github.com/vgteam/vg_wdl/raw/master/workflows/vg_multi_map_call.wdl)
- [parameter file](https://github.com/vgteam/vg_wdl/raw/master/params/vg_multi_map_call.inputs_tiny.http_url.json)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/vg-pipeline-workingexample:master?tab=info)

Workflow for processing pedigree datasets:
- [workflow file](https://github.com/vgteam/vg_wdl/raw/master/workflows/vg_trio_multi_map_call.wdl)
- [parameter file](https://github.com/vgteam/vg_wdl/raw/master/params/vg_trio_multi_map_call.inputs_tiny.http_url.json)

Workflow for mapping and calling structural variants in a single sample:
- [workflow file](https://github.com/vgteam/vg_wdl/raw/svpack/workflows/vg_map_call_sv.wdl)
- [parameter file](https://github.com/vgteam/vg_wdl/raw/svpack/params/vg_map_call_sv_test.inputs.json)
- [Dockstore page](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/vg_map_call_sv:svpack?tab=info)

Workflow for mapping short reads with vg Giraffe and calling short variants with DeepVariant:
- [workflow file](workflows/giraffe_and_deepvariant.wdl)

## Contributing, Help, Bugs and Requests
Please open an Issue on [github](https://github.com/vgteam/vg_wdl) for help, bug reports, or feature requests.
When doing so, please remember that vg\_wdl is open-source software made by a community of developers. Please
be considerate and support a positive environment.
