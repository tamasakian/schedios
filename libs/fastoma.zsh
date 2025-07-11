#!/usr/bin/env zsh

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