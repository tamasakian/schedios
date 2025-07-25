#!/usr/bin/env zsh

function detect_collinearity_within_gene_families() {
    # This function runs downstream analyses for MCScanX collinearity files.
    # It processes each collinearity file in the input directory and generates output files.
    # Dependencies: perl, MCScanX

    # === Usage ===
    if (( $# != 2 )); then
        echo "Usage: run_mcscanx_da <taskname> <family>" >&2
        exit 1
    fi

    # === Arguments ===
    local taskname="$1"
    local family="$2"
    local taskdir="${TASKS}/${taskname}"

    # === Make directories ===
    local outdir="${taskdir}/output/$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "$outdir"

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
    for input_file in "${TASKFILE}/${taskname}/input/collinearity/"*.collinearity; do
        local base_name="${input_file##*/}"
        local output_file="${outdir}/${base_name%.collinearity}.txt"
        perl ${HOME}/bin/MCScanX-1.0.0/downstream_analyses/detect_collinearity_within_gene_families.pl -i "$family_file" -d "$input_file" -o "$output_file"
    done
}