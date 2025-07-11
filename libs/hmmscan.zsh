#!/usr/bin/env zsh

function hmmscan_pfam_domain() {
    # This function runs hmmscan on a given task directory to scan for Pfam domains.
    # dependencies: hmmer

    # === Defaults ===
    local threads=4
    local cutoffs="ga,nc,tc"

    # === Usage ===
    if (( $# < 1 )); then
        echo "Usage: scan_pfam_domain <taskname> [--threads <threads>] [--cutoffs <comma_separated_list_of_cutoffs>]" >&2
        exit 1
    fi

    # === Arguments ===
    local taskname="$1"
    shift

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
    local hmm="${DATA}/Pfam/Pfam-A.hmm"

    # === Make directories ===
    local outdir="${taskdir}/output_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "$outdir"

    # === HMMSCAN Execution ===
    IFS=',' read -r -a cutoff_arr <<< "$cutoffs"
    for cut in "${cutoff_arr[@]}"; do
        hmmscan \
            -o "${outdir}/output/hmmscan_${cut}.txt" \
            --domtblout "${taskdir}/output/domtblout_${cut}.txt" \
            --cpu "$threads" \
            --cut_${cut} \
            "$hmm" \
            "$fasta"
    done
}
