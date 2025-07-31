#!/usr/bin/env zsh

function input4orthofinder(){
    # This function prepares input for OrthoFinder by selecting longest isoforms and converting sequence IDs.
    # Dependencies: bithon, fasp

    # === Arguments ===
    if (( $# < 2 )); then
        echo "Usage: input4orthofinder <sp_name1> <sp_name2> ..." >&2
        exit 1
    fi
    names=("${@}")

    # === Make directories ===
    taskdir="${TASKS}/orthofinder_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/input"

    # === Proteomes ===
    for name in "${names[@]}"; do
        sp_name="${name// /_}"
        pep="${DATA}/${sp_name}/${sp_name}.pep.all.fasta"
        cds="${DATA}/${sp_name}/${sp_name}.cds.all.fasta"
        tmpdir="${taskdir}/tmp_${sp_name}"
        mkdir -p "$tmpdir"
        cp "$pep" "${tmpdir}/pep.fasta"
        cp "$cds" "${tmpdir}/cds.fasta"

        # === Longest isoforms selection with bithon ===
        bithon gls -i "$tmpdir" -o "$tmpdir"

        python3 -m fasp prefix_to_sequence_ids \
            "${tmpdir}/longest.pep.fasta" \
            "${taskdir}/input/${sp_name}.fasta" \
            "$sp_name"

        rm -r "$tmpdir"

    done
}


function run_orthofinder() {
    # This function runs OrthoFinder on a given task directory.
    # Dependencies: orthofinder

    # === Defaults ===
    threads=4

    # === Arguments ===
    if (( $# < 1 )); then
        echo "Usage: run_orthofinder <taskname> [--threads <threads>]" >&2
        exit 1
    fi

    taskname="$1"
    shift

    while (( $# > 0 )); do
        case "$1" in
            --threads)
                threads="$2"
                shift 2
                ;;
            *)
                echo "[Error] Unknown option: $1" >&2
                exit 1
                ;;
        esac
    done

    # === Make directories ===
    taskdir="${TASKS}/${taskname}"

    # === Orthofinder execution ===
    orthofinder -f "${taskdir}/input" -t "$threads"
}