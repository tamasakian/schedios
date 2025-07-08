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
