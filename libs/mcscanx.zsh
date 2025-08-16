#!/usr/bin/env zsh

function input4mcscanx() {
    # This function prepares input for MCScanX by selecting longest isoforms,
    # converting sequence IDs, and preparing BED files for pairwise comparisons.
    # Dependencies: bithon, fasp, g2bp, diamond, MCScanX

    # === Arguments ===
    local threads=4
    if [[ "$1" == "--threads" ]]; then
        threads="$2"
        shift 2
    fi

    if (( $# < 2 )); then
        echo "Usage: input4mcscanx [--threads N] <sp1> <sp2> ..." >&2
        exit 1
    fi
    local names=("$@")

    # === Prepare task directory ===
    local taskdir="${TASKS}/mcscanx_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/fasta" "${taskdir}/bed" "${taskdir}/pairs" "${taskdir}/collinearity"

    # === Prepare files for each species ===
    for name in "${names[@]}"; do
        local sp_name="${name// /_}"
        local tmpdir="${taskdir}/tmp_${sp_name}"
        mkdir -p "$tmpdir"

        local pep="${DATA}/${sp_name}/${sp_name}.pep.all.fasta"
        local cds="${DATA}/${sp_name}/${sp_name}.cds.all.fasta"
        local gff="${DATA}/${sp_name}/${sp_name}.genome.gff"

        cp "$pep" "${tmpdir}/pep.fasta"
        cp "$cds" "${tmpdir}/cds.fasta"
        cp "$gff" "${tmpdir}/genome.gff"

        bithon gls -i "$tmpdir" -o "$tmpdir"

        python3 -m g2bp fasta4mcscanx \
            "${tmpdir}/genome.gff" \
            "${tmpdir}/genome.bed" \
            "${tmpdir}/longest.pep.fasta"

        python3 -m fasp prefix_to_sequence_ids \
            "${tmpdir}/longest.pep.fasta" \
            "${taskdir}/fasta/${sp_name}.fasta" \
            "$sp_name"

        awk -v sp="${sp_name}_" 'BEGIN{OFS=FS="\t"} {$2=sp $2; print}' "${tmpdir}/genome.bed" > "${taskdir}/bed/${sp_name}.bed"

        rm -r "$tmpdir"
    done

    # === Make diamond DB for all species once ===
    for name in "${names[@]}"; do
        local sp_name="${name// /_}"
        diamond makedb --in "${taskdir}/fasta/${sp_name}.fasta" \
            -d "${taskdir}/fasta/${sp_name}.dmnd"
    done

    # === Pairwise comparisons ===
    for ((i=1; i<=${#names[@]}; i++)); do
        for ((j=i+1; j<=${#names[@]}; j++)); do
            local name1="${names[i]}"
            local name2="${names[j]}"
            local sp_name1="${name1// /_}"
            local sp_name2="${name2// /_}"
            local pairdir="${taskdir}/pairs/${sp_name1}__${sp_name2}"
            mkdir -p "$pairdir"

            cat "${taskdir}/bed/${sp_name1}.bed" \
                "${taskdir}/bed/${sp_name2}.bed" \
                > "${pairdir}/${sp_name1}__${sp_name2}.gff"

            diamond blastp \
                -q "${taskdir}/fasta/${sp_name1}.fasta" \
                -d "${taskdir}/fasta/${sp_name2}.dmnd" \
                -o "${pairdir}/${sp_name1}__${sp_name2}.blast" \
                --threads "$threads" --outfmt 6
        done
    done
}

run_mcscanx() {
    # This function runs MCScanX on a given task directory.
    # input4mcscanx should be called before this function to prepare the input.
    # Dependencies: MCScanX

    # === Defaults ===
    local threads=4

    # === Arguments ===
    if (( $# < 1 )); then
        echo "Usage: run_mcscanx <taskname>" >&2
        exit 1
    fi

    local taskname="$1"
    local taskdir="${TASKS}/${taskname}"

    if [[ ! -d "${taskdir}/pairs" ]]; then
        echo "[Error] Input pairs directory not found: ${taskdir}/pairs" >&2
        exit 1
    fi

    if [[ ! -d "${taskdir}/collinearity" ]]; then
        mkdir -p "${taskdir}/collinearity"
    fi

    # === Run MCScanX for each pair ===
    for pairdir in "${taskdir}/pairs/"*; do
        if [[ -d "$pairdir" ]]; then
            local pairname="${pairdir:t}"
            MCScanX "${pairdir}/${pairname}"
            cp "${pairdir}/${pairname}.collinearity" "${taskdir}/collinearity/${pairname}.collinearity"
        fi
    done
}


function input4mcscanx_group() {
    # Prepare input for MCScanX with focus × other comparisons only
    # Dependencies: bithon, fasp, g2bp, diamond, MCScanX

    # === Arguments ===
    local threads=4
    if [[ "$1" == "--threads" ]]; then
        threads="$2"
        shift 2
    fi

    if [[ "$1" != "--group" ]]; then
        echo "Usage: input4mcscanx_group [--threads N] --group <sp1> ... --all <spA> ..." >&2
        return 1
    fi
    shift

    local group1=()
    while [[ "$1" != "--all" && $# -gt 0 ]]; do
        group1+=("$1")
        shift
    done

    if [[ "$1" != "--all" ]]; then
        echo "[Error] Missing --all argument" >&2
        return 1
    fi
    shift

    local all=("$@")

    # === Checks ===
    if (( ${#group1[@]} < 1 )); then
        echo "[Error] At least one group1 species required" >&2
        return 1
    fi
    if (( ${#all[@]} < 2 )); then
        echo "[Error] Need at least two total species" >&2
        return 1
    fi

    # === Determine other species ===
    local group2=()
    for name in "${all[@]}"; do
        local found=0
        for name1 in "${group1[@]}"; do
            [[ "$name" == "$name1" ]] && found=1 && break
        done
        (( found == 0 )) && group2+=("$name")
    done

    # === Prepare task directory ===
    local taskdir="${TASKS}/mcscanx_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/fasta" "${taskdir}/bed" "${taskdir}/pairs" "${taskdir}/collinearity"

    # === Prepare files for each species ===
    for name in "${all[@]}"; do
        local sp_name="${name// /_}"
        local tmpdir="${taskdir}/tmp_${sp_name}"
        mkdir -p "$tmpdir"

        local pep="${DATA}/${sp_name}/${sp_name}.pep.all.fasta"
        local cds="${DATA}/${sp_name}/${sp_name}.cds.all.fasta"
        local gff="${DATA}/${sp_name}/${sp_name}.genome.gff"

        cp "$pep" "${tmpdir}/pep.fasta"
        cp "$cds" "${tmpdir}/cds.fasta"
        cp "$gff" "${tmpdir}/genome.gff"

        bithon gls -i "$tmpdir" -o "$tmpdir"

        python3 -m g2bp fasta4mcscanx \
            "${tmpdir}/genome.gff" \
            "${tmpdir}/genome.bed" \
            "${tmpdir}/longest.pep.fasta"

        python3 -m fasp prefix_to_sequence_ids \
            "${tmpdir}/longest.pep.fasta" \
            "${taskdir}/fasta/${sp_name}.fasta" \
            "$sp_name"

        awk -v sp="${sp_name}_" 'BEGIN{OFS=FS="\t"} {$2=sp $2; print}' "${tmpdir}/genome.bed" > "${taskdir}/bed/${sp_name}.bed"

        rm -r "$tmpdir"
    done

    # === Make diamond DB for group1 species once ===
    for name in "${group1[@]}"; do
        local sp_name="${name// /_}"
        diamond makedb --in "${taskdir}/fasta/${sp_name}.fasta" \
            -d "${taskdir}/fasta/${sp_name}.dmnd"
    done

    # === Pairwise comparisons: group1 × group2 ===
    for name1 in "${group1[@]}"; do
        for name2 in "${group2[@]}"; do
            local sp_name1="${name1// /_}"
            local sp_name2="${name2// /_}"
            local pairdir="${taskdir}/pairs/${sp_name1}__${sp_name2}"
            mkdir -p "$pairdir"

            cat "${taskdir}/bed/${sp_name1}.bed" \
                "${taskdir}/bed/${sp_name2}.bed" \
                > "${pairdir}/${sp_name1}__${sp_name2}.gff"

            diamond blastp \
                -q "${taskdir}/fasta/${sp_name2}.fasta" \
                -d "${taskdir}/fasta/${sp_name1}.dmnd" \
                -o "${pairdir}/${sp_name1}__${sp_name2}.blast" \
                --threads "$threads" --outfmt 6
        done
    done
}



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