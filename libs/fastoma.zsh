#!/usr/bin/env zsh

function input4fastoma(){
    # This function prepares input for FastOMA by selecting longest isoforms and converting sequence IDs.
    # Dependencies: bithon, biotp, fasp, ete3

    # === Arguments ===
    if (( $# < 2 )); then
        echo "Usage: input4fastoma <sp_name1> <sp_name2> ..." >&2
        exit 1
    fi
    names=("${@}")

    # === Make directories ===
    taskdir="${TASKS}/fastoma_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "${taskdir}/proteomes"

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
            "${taskdir}/proteomes/${sp_name}.fasta" \
            "$sp_name"

        rm -r "$tmpdir"

    done

    # === Species tree ===
    print -l "${names[@]}" > "${taskdir}/names.txt"
    taxonkit name2taxid < "${taskdir}/names.txt" | cut -f2 > "${taskdir}/list_taxanomic_id.txt"
    cat "${taskdir}/list_taxanomic_id.txt" | ete3 ncbiquery --tree > "${taskdir}/sp_tree.nwk"
    python3 -m biotp clean_tip_label "${taskdir}/sp_tree.nwk" "${taskdir}/species_tree.nwk"
}

function run_fastoma() {
    # This function runs FastOMA on a given task directory.
    # Dependencies: fastoma

    # === Defaults ===
    local threads=4
    local memory=8.GB
    local db="LUCA"

    # === Usage ===
    if (( $# < 1 )); then
        echo "Usage: run_fastoma <taskname> [--threads <threads>] [--memory <memory>] [--db <db>]" >&2
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
            --memory)
                memory="$2"
                shift 2
                ;;
            --db)
                db="$2"
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
    local indir="${taskdir}/input/"

    # === Make directories ===
    local outdir="${taskdir}/output_$(date +"%Y-%m-%d-%H-%M-%S")"
    mkdir -p "$outdir"

    # === FastOMA Execution ===
    nextflow run dessimozlab/FastOMA \
        -r v0.3.5 \
        -profile singularity \
        --input_folder "$indir" \
        --output_folder "$outdir" \
        --omamer_db "${DATA}/OMAmer/${db}.h5" \
        -process.cpus=${threads} \
        -process.memory=${memory}
}