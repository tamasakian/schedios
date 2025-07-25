#!/usr/bin/env zsh

function detect_collinearity_within_gene_families() {
    # This function runs downstream analyses for MCScanX collinearity files.
    # It processes each collinearity file in the input directory and generates output files.
    # Dependencies: perl, MCScanX

    # === Usage ===
    if (( $# != 2 )); then
        echo "Usage: detect_collinearity_within_gene_families <taskname> <family>" >&2
        exit 1
    fi

    # === Arguments ===
    local taskname="$1"
    local family="$2"
    local taskdir="${TASKS}/${taskname}"

    # === Make directories ===
    local outdir="${taskdir}/output/$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/output" "$outdir" "${outdir}/synteny" "${outdir}/long"

    # === Check input directory ===
    if [[ ! -d "${taskdir}/input" ]]; then
        echo "[Error] Input directory not found: ${taskdir}/input" >&2
        exit 1
    fi

    # === Check collinearity files ===
    if [[ ! -d "${taskdir}/input/collinearity/" ]]; then
        echo "[Error] Input collinearity directory not found: ${taskdir}/input/collinearity/" >&2
        exit 1
    fi

    # === Process each collinearity file ===
    local family_file="${taskdir}/input/${family}"

    if [[ ! -f "$family_file" ]]; then
        echo "[Error] Family file not found: $family_file" >&2
        exit 1
    fi

    local sog_file="${outdir}/sog.txt"
    touch "$sog_file"
    for input_file in "${taskdir}/input/collinearity/"*.collinearity; do
        local base_name="${input_file##*/}"
        local output_file="${outdir}/synteny/${base_name%.collinearity}.txt"
        # === Run the Perl script to detect collinearity within gene families ===
        echo "[INFO] Processing collinearity file: $input_file"
        perl ${HOME}/bin/MCScanX-1.0.0/downstream_analyses/detect_collinearity_within_gene_families.pl -i "$family_file" -d "$input_file" -o "$output_file"

        # === Convert output to long format and append to SOG file ===
        echo "[INFO] Converting output to long format: $output_file"
        python3 -m tidyg pivot_longer "$output_file" "${outdir}/long/${base_name%.collinearity}.txt"
        cat "${outdir}/long/${base_name%.collinearity}.txt" >> "$sog_file"
    done
}