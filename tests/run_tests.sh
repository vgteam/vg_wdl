#!/bin/bash
set -ex -o pipefail

ORIG_CWD="$(pwd)"
cd "$(dirname $0)/ABOlocus"
#miniwdl check --no-quant-check vg_ABOlocus_test_small.wdl
#miniwdl cromwell --no-quant-check -d "$ORIG_CWD/vg_wdl_tests_ABOlocus_small" \
#    vg_ABOlocus_test_small.wdl \
#    ABOlocus_fa_gz=ABOlocus.fa.gz ABOlocus_small_vcf_gz=ABOlocus_small.vcf.gz reads_bam=HG01308.ABOlocus.bam
#miniwdl check --no-quant-check vg_ABOlocus_test_SV.wdl
#miniwdl cromwell --no-quant-check -d "$ORIG_CWD/vg_wdl_tests_ABOlocus_SV" \
#    vg_ABOlocus_test_SV.wdl \
#    ABOlocus_fa_gz=ABOlocus.fa.gz ABOlocus_SV_vcf_gz=ABOlocus_SV.vcf.gz reads_bam=HG01308.ABOlocus.bam
miniwdl check --no-quant-check vg_ABOlocus_test_pedigree.wdl
miniwdl cromwell --no-quant-check -d "$ORIG_CWD/vg_wdl_tests_ABOlocus_pedigree" \
    vg_ABOlocus_test_pedigree.wdl \
    ABOlocus_fa_gz=ABOlocus.fa.gz ABOlocus_small_vcf_gz=ABOlocus_small.vcf.gz \
    maternal_reads_bam=HG004.hs37d5.2x250.abo.bam \
    paternal_reads_bam=HG003.hs37d5.2x250.abo.bam \
    child_reads_bam=HG002.hs37d5.2x250.abo.bam \
    ref_file=ABOlocus.fa \
    ref_index_file=ABOlocus.fa.fai \
    ref_dict_file=ABOlocus.dict \
    ref_file_gz=ABOlocus.fa.gz \
    ped_file=HG002.ped
