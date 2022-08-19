version 1.0

import "../tasks/gam_gaf_utils.wdl" as gautils

workflow sortGraphAlignedReads {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Sort a GAM file or a GAF file. Inputs are either: 1) a GAM file, 2)a GAF file and the XG index."
    }
    parameter_meta {
        GAM_FILE: "GAM file to sort."
        GAF_FILE: "GAF file to sort. The xg index must also be passed in XG_FILE"
        XG_FILE: "the xg index of the graph (necessary when the input if a GAF file)"
    }
    
    input {
        File? GAM_FILE
        File? GAF_FILE
        File? XG_FILE
        String? SAMPLE_NAME
    }

    if(defined(GAM_FILE)) {
	    call gautils.mergeGAMandSort {
            input:
            in_gam_files=select_all([GAM_FILE]),
            in_sample_name=SAMPLE_NAME
	    }
    }

    if(defined(GAF_FILE) && defined(XG_FILE)) {
	    call gautils.mergeGAFandSort {
            input:
            in_gaf_files=select_all([GAF_FILE]),
            in_xg_file=select_first([XG_FILE]),
            in_sample_name=SAMPLE_NAME
	    }
    }

    File out_sorted_gam = select_first([mergeGAMandSort.gam, mergeGAFandSort.gam])
    File out_sorted_gam_index = select_first([mergeGAMandSort.gam_index, mergeGAFandSort.gam_index])
    
    output {
        File sorted_gam = out_sorted_gam
        File sorted_gam_index = out_sorted_gam_index
        }   
    }
    
