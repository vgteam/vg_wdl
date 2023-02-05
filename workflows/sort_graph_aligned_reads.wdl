version 1.0

import "../tasks/gam_gaf_utils.wdl" as gautils

workflow sortGraphAlignedReads {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Convert a GAF file make a sorted/indexed GAM. More information at [https://github.com/vgteam/vg_wdl/tree/gbz#gaf-to-sorted-gam-workflow](https://github.com/vgteam/vg_wdl/tree/gbz#gaf-to-sorted-gam-workflow)."
    }

    parameter_meta {
        GAF_FILE: "GAF file to convert and sort. "
        GBZ_FILE: "the GBZ index of the graph"
        SAMPLE_NAME: "(Optional) a sample name"
    }
    
    input {
        File GAF_FILE
        File GBZ_FILE
        String? SAMPLE_NAME
    }

	call gautils.mergeGAFandSort {
        input:
        in_gaf_file=GAF_FILE,
        in_gbz_file=GBZ_FILE,
        in_sample_name=SAMPLE_NAME
	}
    
    output {
        File sorted_gam = mergeGAFandSort.gam
        File sorted_gam_index = mergeGAFandSort.gam_index
        }   
    }
    
