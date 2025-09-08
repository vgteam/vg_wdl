#!/usr/bin/env bash
set -euo pipefail

declare -A WDL_MAP

# Manually define JSON → WDL mappings
WDL_MAP["giraffe.json"]="giraffe.wdl"
WDL_MAP["giraffe.singleended.cram.json"]="giraffe.wdl"
WDL_MAP["giraffe_and_deepvariant.json"]="giraffe_and_deepvariant.wdl"
WDL_MAP["giraffe_and_deepvariant_cram.json"]="giraffe_and_deepvariant.wdl"
WDL_MAP["giraffe_and_deepvariant_gaf.json"]="giraffe_and_deepvariant_fromGAF.wdl"
WDL_MAP["giraffe_and_deepvariant_gaf_single_end.json"]="giraffe_and_deepvariant_fromGAF.wdl"
WDL_MAP["giraffe_and_deepvariant_single_end.json"]="giraffe_and_deepvariant.wdl"
WDL_MAP["giraffe_and_haplotype_sampling.json"]="giraffe.wdl"
WDL_MAP["haplotype_sampling.json"]="haplotype_sampling.wdl"
#WDL_MAP["vg_multi_map_call.inputs_21_22.test.json"]="vg_multi_map_call.wdl"
#WDL_MAP["vg_trio_deepvariant_map_call.inputs_tiny.gs_url.json"]="vg_trio_deepvariant_map_call.wdl"
# Add the rest of your mappings here...

# Track failures
FAILURES=()

echo "Running mapped WDL workflows..."

for json_file in "${!WDL_MAP[@]}"; do
    input_path="params/${json_file}"
    wdl_path="workflows/${WDL_MAP[$json_file]}"

    echo "---------------------------------------------"
    echo "Running ${wdl_path} with ${input_path}"

    if [[ ! -f "$input_path" ]]; then
        echo "❌ Missing input JSON: $input_path"
        FAILURES+=("$json_file (missing input)")
        continue
    fi

    if [[ ! -f "$wdl_path" ]]; then
        echo "❌ Missing WDL file: $wdl_path"
        FAILURES+=("$json_file (missing WDL)")
        continue
    fi

    if miniwdl run "$wdl_path" --input "$input_path"; then
        echo "✅ Success: $wdl_path with $json_file"
    else
        echo "❌ Failure: $wdl_path with $json_file"
        FAILURES+=("$json_file")
    fi
done

echo "============================================="
if [ ${#FAILURES[@]} -eq 0 ]; then
    echo "🎉 All mapped workflows ran successfully!"
    exit 0
else
    echo "❌ Some workflows failed:"
    for f in "${FAILURES[@]}"; do
        echo " - $f"
    done
    exit 1
fi
