#!/usr/bin/env zsh

function forge_hmmscan_input() {
    # This function prepares input for hmmscan by creating a task directory
    # and generating a combined FASTA file from multiple species.
    # Dependencies: bithon, fasp

    # === Usage ===
    if (( $# < 2 )); then
        echo "Usage: forge_pfam_input <sp_name1> <sp_name2> ..." >&2
        exit 1
    fi

    # === Arguments ===
    local sp_names=("${@}")

    # === Make directories ===
    local taskdir="${TASKS}/hmmscan_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/tmp" "${taskdir}/fasta" "${taskdir}/input"

    for sp in "${sp_names[@]}"; do
        local sp_fs="${sp// /_}"
        local gn="${sp_fs%%_*}"
        local pep="${DATA}/${gn}/${sp_fs}.pep.all.fasta"
        local cds="${DATA}/${gn}/${sp_fs}.cds.all.fasta"

        mkdir -p "${taskdir}/tmp/${sp_fs}"
        cp "${pep}" "${taskdir}/tmp/${sp_fs}/pep.fasta"
        cp "${cds}" "${taskdir}/tmp/${sp_fs}/cds.fasta"

        # === Longest isoforms selection with bithon ===
        echo "[INFO] Processing species: ${sp} (Pep: ${pep}, CDS: ${cds})"
        bithon gls -i "${taskdir}/tmp/${sp_fs}" -o "${taskdir}/tmp/${sp_fs}"
        cp "${taskdir}/tmp/${sp_fs}/longest.pep.fasta" "${taskdir}/fasta/${sp_fs}.fasta"
        rm -r "${taskdir}/tmp/${sp_fs}"

        # === Prefix to sequence IDs conversion ===
        echo "[INFO] Converting sequence IDs for species: ${sp}"
        python3 -m fasp prefix_to_sequence_ids \
            "${taskdir}/fasta/${sp_fs}.fasta" \
            "${taskdir}/fasta/${sp_fs}_prefixed.fasta" \
            "${sp_fs}"
    done

    # === Combine all FASTA files into one ===
    echo "[FINISHED] HMMscan input prepared in: ${taskdir}"
    cat "${taskdir}/fasta/"*_prefixed.fasta > "${taskdir}/input/all_species.fasta"
}

function forge_fastoma_input(){
    # This function prepares input for FastOMA by creating a task directory.
    # Dependencies: bithon, biotp, fasp, ete3

    # === Usage ===
    if (( $# < 2 )); then
        echo "Usage: forge_fastoma_input <sp_name1> <sp_name2> ..." >&2
        exit 1
    fi

    # === Arguments ===
    local sp_names=("${@}")

    # === Make directories ===
    local taskdir="${TASKS}/fastoma_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/tmp" "${taskdir}/fasta" "${taskdir}/input" "${taskdir}/input/proteome"

    for sp in "${sp_names[@]}"; do
        local sp_fs="${sp// /_}"
        local gn="${sp_fs%%_*}"
        local pep="${DATA}/${gn}/${sp_fs}.pep.all.fasta"
        local cds="${DATA}/${gn}/${sp_fs}.cds.all.fasta"

        mkdir -p "${taskdir}/tmp/${sp_fs}"
        cp "${pep}" "${taskdir}/tmp/${sp_fs}/pep.fasta"
        cp "${cds}" "${taskdir}/tmp/${sp_fs}/cds.fasta"

        # === Longest isoforms selection with bithon ===
        echo "[INFO] Processing species: ${sp} (Pep: ${pep}, CDS: ${cds})"
        bithon gls -i "${taskdir}/tmp/${sp_fs}" -o "${taskdir}/tmp/${sp_fs}"
        cp "${taskdir}/tmp/${sp_fs}/longest.pep.fasta" "${taskdir}/fasta/${sp_fs}.fasta"
        rm -r "${taskdir}/tmp/${sp_fs}"

        # === Prefix to sequence IDs conversion ===
        echo "[INFO] Converting sequence IDs for species: ${sp}"
        python3 -m fasp prefix_to_sequence_ids \
            "${taskdir}/fasta/${sp_fs}.fasta" \
            "${taskdir}/input/proteome/${sp_fs}.fasta" \
            "${sp_fs}"
    done

    # === Species tree generation ===
    echo "[INFO] Generating species tree"
    local ete_args=()
    for sp in "${sp_names[@]}"; do
        ete_args+=("'${sp}'")
    done
    local ete_arg=$(IFS=" "; echo "${ete_args[*]}")
    newick_str=$(ete3 ncbiquery --search ${ete_arg} --tree)
    echo "$newick_str" | python3 -m biotp clean_leaf_names > "${taskdir}/input/species_tree.nwk"
    echo "[FINISHED] FastOMA input prepared in: ${taskdir}"
}

function forge_mcscanx() {
    # This function prepares input for MCScanX by creating a task directory
    # and generating pairwise comparisons of species.
    # Dependencies: bithon, fasp, diamond, MCScanX

    # === Defaults ===
    local threads=4
    if [[ "$1" == "--threads" ]]; then
        threads="$2"
        shift 2
    fi

    # === Usage ===
    if (( $# < 1 )); then
        echo "Usage: forge_mcscanx_input [--threads <threads> <sp_name1> <sp_name2> ..." >&2
        exit 1
    fi

    # === Arguments ===
    local sp_names=("${@}")
    local taskdir="${TASKS}/mcscanx_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/tmp" "${taskdir}/fasta" "${taskdir}/fasta_prefix" "${taskdir}/gff" "${taskdir}/bed" "${taskdir}/bed_prefix" "${taskdir}/input"
    mkdir -p "${taskdir}/input/pairs" "${taskdir}/input/collinearity"

    echo "[INFO] Preparing MCScanX input for species: ${sp_names[*]}"
    if [[ ${#sp_names[@]} -lt 2 ]]; then
        echo "[Error] At least two species are required for MCScanX." >&2
        exit 1
    fi

    # === Prepare FASTA, GFF, and BED files for each species ===
    for sp in "${sp_names[@]}"; do
        echo "[INFO] Processing species: ${sp}"
        local sp_fs="${sp// /_}"
        local gn="${sp_fs%%_*}"

        local pep="${DATA}/${gn}/${sp_fs}.pep.all.fasta"
        local cds="${DATA}/${gn}/${sp_fs}.cds.all.fasta"
        local gff="${DATA}/${gn}/${sp_fs}.genome.gff"

        mkdir -p "${taskdir}/tmp/${sp_fs}"
        cp "$pep" "${taskdir}/tmp/${sp_fs}/pep.fasta"
        cp "$cds" "${taskdir}/tmp/${sp_fs}/cds.fasta"
        cp "$gff" "${taskdir}/gff/${sp_fs}.gff"
        bithon gls -i "${taskdir}/tmp/${sp_fs}" -o "${taskdir}/tmp/${sp_fs}"
        cp "${taskdir}/tmp/${sp_fs}/longest.pep.fasta" "${taskdir}/fasta/${sp_fs}.fasta"
        rm -r "${taskdir}/tmp/${sp_fs}"

        python3 -m g2bp fasta4mcscanx "${taskdir}/gff/${sp_fs}.gff" "${taskdir}/bed/${sp_fs}.bed" "${taskdir}/fasta/${sp_fs}.fasta"
        python3 -m fasp prefix_to_sequence_ids "${taskdir}/fasta/${sp_fs}.fasta" "${taskdir}/fasta_prefix/${sp_fs}.fasta" "$sp_fs"
        awk -v sp="${sp_fs}_" 'BEGIN{OFS=FS="\t"} {$2=sp $2; print}' "${taskdir}/bed/${sp_fs}.bed" > "${taskdir}/bed_prefix/${sp_fs}.bed"

    done

    # === Create pairs and run MCScanX ===
    for ((i=1; i<=${#sp_names[@]}-1; i++)); do
        for ((j=i+1; j<=${#sp_names[@]}; j++)); do
            local sp1="${sp_names[i]}"
            local sp2="${sp_names[j]}"
            local sp1_fs="${sp1// /_}"
            local sp2_fs="${sp2// /_}"
            local gn1="${sp1_fs%%_*}"
            local gn2="${sp2_fs%%_*}"
            local pairdir="${taskdir}/input/pairs/${sp1_fs}__${sp2_fs}"
            mkdir -p "$pairdir"

            echo "[INFO] Processing pair: ${sp1} vs ${sp2}"
            cp "${taskdir}/fasta_prefix/${sp1_fs}.fasta" "${pairdir}/${sp1_fs}.fasta"
            cp "${taskdir}/fasta_prefix/${sp2_fs}.fasta" "${pairdir}/${sp2_fs}.fasta"
            cat "${taskdir}/bed_prefix/${sp1_fs}.bed" "${taskdir}/bed_prefix/${sp2_fs}.bed" > "${pairdir}/${sp1_fs}__${sp2_fs}.gff"
            diamond makedb --in "${pairdir}/${sp2_fs}.fasta" -d "${pairdir}/${sp2_fs}.dmnd"
            diamond blastp \
                -q "${pairdir}/${sp1_fs}.fasta" \
                -d "${pairdir}/${sp2_fs}.dmnd" \
                -o "${pairdir}/${sp1_fs}__${sp2_fs}.blast" \
                --threads "$threads" \
                --outfmt 6
            MCScanX "${pairdir}/${sp1_fs}__${sp2_fs}"
            cp "${pairdir}/${sp1_fs}__${sp2_fs}.collinearity" "${taskdir}/input/collinearity/${sp1_fs}__${sp2_fs}.collinearity"
        done
    done
}

