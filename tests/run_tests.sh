#!/bin/bash
set -ex -o pipefail

ORIG_CWD="$(pwd)"
cd "$(dirname $0)/ABOlocus"
miniwdl check --no-quant-check vg_ABOlocus_test_small.wdl
miniwdl cromwell --no-quant-check -d "$ORIG_CWD/vg_wdl_tests_ABOlocus_small" vg_ABOlocus_test_small.wdl \
    ABOlocus_fa_gz=ABOlocus.fa.gz ABOlocus_small_vcf_gz=ABOlocus_small.vcf.gz reads_bam=HG01308.ABOlocus.bam
miniwdl check --no-quant-check vg_ABOlocus_test_SV.wdl
miniwdl cromwell --no-quant-check -d "$ORIG_CWD/vg_wdl_tests_ABOlocus_SV" vg_ABOlocus_test_SV.wdl \
    ABOlocus_fa_gz=ABOlocus.fa.gz ABOlocus_SV_vcf_gz=ABOlocus_SV.vcf.gz reads_bam=HG01308.ABOlocus.bam
