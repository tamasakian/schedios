#!/usr/bin/env zsh

function hmmsearch_pfam_domain() {
    # This function runs hmmsearch on a given task directory to search for Pfam domains.
    # dependencies: hmmer

    # === Defaults ===
    local threads=4
    local cutoffs="ga,nc,tc"

    # === Usage ===
    if (( $# < 1 )); then
        echo "Usage: hmmsearch_pfam_domain <taskname> <domain> [--threads <threads>] [--cutoffs <comma_separated_list_of_cutoffs>]" >&2
        exit 1
    fi

    # === Arguments ===
    local taskname="$1"
    local domain="$2"
    shift 2

    while (( $# > 0 )); do
        case "$1" in
            --threads)
                threads="$2"
                shift 2
                ;;
            --cutoffs)
                cutoffs="$2"
                shift 2
                ;;
            *)
                echo "[Error] Unknown option: $1" >&2
                exit 1
                ;;
        esac
    done

    # === Set paths ===
    local taskdir="${TASKS}/${taskname}"
    local fasta="${taskdir}/input/all_species.fasta"
    local hmm="${DATA}/Pfam/${domain}.hmm"

    # === Make directories ===
    local outdir="${taskdir}/output_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "$outdir"

    # === HMMSEARCH Execution ===
    cutoff_arr=("${(s:,:)cutoffs}")
    for cut in "${cutoff_arr[@]}"; do
        hmmsearch \
            -o "${outdir}/hmmsearch_${cut}.txt" \
            --domtblout "${outdir}/domtblout_${cut}.txt" \
            --cpu "$threads" \
            --cut_${cut} \
            "$hmm" \
            "$fasta"
    done
}