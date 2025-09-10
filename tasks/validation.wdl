version 1.0

# validation.wdl: Tasks to check metadata files like sequence dictionaries and path lists to make sure they obey workflow requirements

# Ensure that all paths in the path list have the provided reference prefix.
task checkPathList {
    input {
        File in_path_list_file
        String in_reference_prefix
    }

    command <<<
        set -eu

        if [[ "~{in_reference_prefix}" == "" ]] ; then
            # No prefix used
            exit 0
        fi
        
        # To match a prefix and also use a non-regex fixed string pattern, we
        # cut each line to length like StackOverflow recommends.
        UNPREFIXED_COUNT="$(cat ~{in_path_list_file} | cut -c 1-$(echo -n "~{in_reference_prefix}" | wc -c) | grep -v -F "~{in_reference_prefix}" | wc -l)"
        if [[ "${UNPREFIXED_COUNT}" != "0" ]] ; then
            echo 1>&2 "Error: path list file ~{in_path_list_file} contains paths that do not start with reference prefix \"~{in_reference_prefix}\"."
            exit 1
        fi
    >>>

    output {
        # No output needed; we fail if we don't like the file
    }

    runtime {
        preemptible: 2
        time: 20
        cpu: 1
        memory: "2 GB"
        disks: "local-disk " + (size(in_path_list_file, "GB") + 1) + " SSD"
        docker: "ubuntu:24.04"
    }
}

# Ensure that no path in the dict file has the reference prefix
task checkDict {
    input {
        File in_dict_file
        String in_reference_prefix
    }

    command <<<
        set -eu

        if [[ "~{in_reference_prefix}" == "" ]] ; then
            # No prefix used
            exit 0
        fi
        
        # To match a prefix and also use a non-regex fixed string pattern, we
        # cut each line to length like StackOverflow recommends.
        PREFIXED_COUNT="$(cat ~{in_dict_file} | grep "^@SQ" | tr '\t' '\n' | grep "^SN:" | cut -f2 -d':' | cut -c 1-$(echo -n "~{in_reference_prefix}" | wc -c) | grep -F "~{in_reference_prefix}" | wc -l)"
        if [[ "${PREFIXED_COUNT}" != "0" ]] ; then
            echo 1>&2 "Error: dict file ~{in_dict_file} contains paths that start with reference prefix \"~{in_reference_prefix}\"."
            exit 1
        fi
    >>>

    output {
        # No output needed; we fail if we don't like the file
    }

    runtime {
        preemptible: 2
        time: 20
        cpu: 1
        memory: "2 GB"
        disks: "local-disk " + (size(in_dict_file, "GB") + 1) + " SSD"
        docker: "ubuntu:24.04"
    }

}
