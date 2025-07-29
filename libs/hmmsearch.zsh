#!/usr/bin/env zsh

function forge_hmmsearch_input() {
    # This function prepares input for HMMsearch by creating a task directory
    # and generating a combined FASTA file from multiple species.
    # Dependencies: bithon, fasp

    # === Arguments ===
    if (( $# < 2 )); then
        echo "Usage: forge_hmmsearch_input <sp_name1> <sp_name2> ..." >&2
        exit 1
    fi
    names=("${@}")

    # === Make directories ===
    taskdir="${TASKS}/hmmsearch_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/fasta" "${taskdir}/input"

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
            "${taskdir}/fasta/${sp_name}.fasta" \
            "$sp_name"

        rm -r "$tmpdir"

    done

    cat "${taskdir}/fasta/"*.fasta > "${taskdir}/input/all_species.fasta"
}

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