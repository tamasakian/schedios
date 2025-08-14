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
        return 1
    fi
    local all=("$@")

    # === Prepare task directory ===
    local taskdir="${TASKS}/mcscanx_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/fasta" "${taskdir}/bed" "${taskdir}/pairs" "${taskdir}/collinearity"

    # === Prepare files for each species ===
    for sp in "${all[@]}"; do
        local sp_name="${sp// /_}"
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

    # === Make diamond DB for focus species once ===
    for sp in "${all[@]}"; do
        local sp_name="${sp// /_}"
        diamond makedb --in "${taskdir}/fasta/${sp_name}.fasta" \
            -d "${taskdir}/fasta/${sp_name}.dmnd"
    done

    # === Pairwise comparisons: focus × other ===
    for ((i=0; i<${#all[@]}; i++)); do
        for ((j=i+1; j<${#all[@]}; j++)); do
            local sp1="${all[i]}"
            local sp2="${all[j]}"
            local sp1_name="${sp1// /_}"
            local sp2_name="${sp2// /_}"
            local pairdir="${taskdir}/pairs/${sp1_name}__${sp2_name}"
            mkdir -p "$pairdir"

            cat "${taskdir}/bed/${sp1_name}.bed" \
                "${taskdir}/bed/${sp2_name}.bed" \
                > "${pairdir}/${sp1_name}__${sp2_name}.gff"

            diamond blastp \
                -q "${taskdir}/fasta/${sp1_name}.fasta" \
                -d "${taskdir}/fasta/${sp2_name}.dmnd" \
                -o "${pairdir}/${sp1_name}__${sp2_name}.blast" \
                --threads "$threads" --outfmt 6

            MCScanX "${pairdir}/${sp1_name}__${sp2_name}"

            cp "${pairdir}/${sp1_name}__${sp2_name}.collinearity" "${taskdir}/collinearity/${sp1_name}__${sp2_name}.collinearity"
        done
    done
}

function input4mcscanx2() {
    # Prepare input for MCScanX with focus × other comparisons only
    # Dependencies: bithon, fasp, g2bp, diamond, MCScanX

    # === Arguments ===
    local threads=4
    if [[ "$1" == "--threads" ]]; then
        threads="$2"
        shift 2
    fi

    if [[ "$1" != "--focus" ]]; then
        echo "Usage: input4mcscanx2 [--threads N] --focus <sp1> ... --all <spA> ..." >&2
        return 1
    fi
    shift

    local focus=()
    while [[ "$1" != "--all" && $# -gt 0 ]]; do
        focus+=("$1")
        shift
    done

    if [[ "$1" != "--all" ]]; then
        echo "[Error] Missing --all argument" >&2
        return 1
    fi
    shift

    local all=("$@")

    # === Checks ===
    if (( ${#focus[@]} < 1 )); then
        echo "[Error] At least one focus species required" >&2
        return 1
    fi
    if (( ${#all[@]} < 2 )); then
        echo "[Error] Need at least two total species" >&2
        return 1
    fi

    # === Determine other species ===
    local others=()
    for sp in "${all[@]}"; do
        local found=0
        for fsp in "${focus[@]}"; do
            [[ "$sp" == "$fsp" ]] && found=1 && break
        done
        (( found == 0 )) && others+=("$sp")
    done

    # === Prepare task directory ===
    local taskdir="${TASKS}/mcscanx_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/fasta" "${taskdir}/bed" "${taskdir}/pairs" "${taskdir}/collinearity"

    # === Prepare files for each species ===
    for sp in "${all[@]}"; do
        local sp_name="${sp// /_}"
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

    # === Make diamond DB for focus species once ===
    for sp in "${focus[@]}"; do
        local sp_name="${sp// /_}"
        diamond makedb --in "${taskdir}/fasta/${sp_name}.fasta" \
            -d "${taskdir}/fasta/${sp_name}.dmnd"
    done

    # === Pairwise comparisons: focus × other ===
    for sp1 in "${focus[@]}"; do
        for sp2 in "${others[@]}"; do
            local sp1_name="${sp1// /_}"
            local sp2_name="${sp2// /_}"
            local pairdir="${taskdir}/pairs/${sp1_name}__${sp2_name}"
            mkdir -p "$pairdir"

            cat "${taskdir}/bed/${sp1_name}.bed" \
                "${taskdir}/bed/${sp2_name}.bed" \
                > "${pairdir}/${sp1_name}__${sp2_name}.gff"

            diamond blastp \
                -q "${taskdir}/fasta/${sp2_name}.fasta" \
                -d "${taskdir}/fasta/${sp1_name}.dmnd" \
                -o "${pairdir}/${sp1_name}__${sp2_name}.blast" \
                --threads "$threads" --outfmt 6

            MCScanX "${pairdir}/${sp1_name}__${sp2_name}"

            cp "${pairdir}/${sp1_name}__${sp2_name}.collinearity" "${taskdir}/collinearity/${sp1_name}__${sp2_name}.collinearity"
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