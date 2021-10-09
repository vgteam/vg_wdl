version 1.0

### vg_deeptrio_calling_workflow.wdl ###
## Author: Charles Markello
## Description: Core VG DEEPTRIO workflow for calling variants on trio datasets.
## Reference: https://github.com/vgteam/vg/wiki

workflow testworkflow {
    input {
        File CHILD_BAM_FILE = "/test/HG002.chr1.indel_realigned.bam"
    }
    String bam_file_basename = basename(CHILD_BAM_FILE)
    String contig_name1 = sub(bam_file_basename, "\\.indel_realigned.bam", "")
    String contig_name2 = sub(sub(bam_file_basename, "\\.indel_realigned.bam", ""), "HG002", "")
    String contig_name3 = sub(contig_name2, "\\.", "")
    String contig_name = sub(sub(sub(bam_file_basename, "\\.indel_realigned.bam", ""), "HG002", ""), "\\.", "")
    if ((contig_name == "chrX")||(contig_name == "X")||(contig_name == "chrY")||(contig_name == "Y")||(contig_name == "chrM")||(contig_name == "MT")) {
        Boolean isNotAutosome = false
    }
    if (!((contig_name == "chrX")||(contig_name == "X")||(contig_name == "chrY")||(contig_name == "Y")||(contig_name == "chrM")||(contig_name == "MT"))) {
        Boolean isAutosome = true
    }
    Boolean isAutosome_result = select_first([isNotAutosome,isAutosome])
    output {
        String contig_name_out = contig_name
        String contig_name_out1 = contig_name1
        String contig_name_out2 = contig_name2
        String contig_name_out3 = contig_name3
        Boolean isAutosome_out = isAutosome_result
    }
}

