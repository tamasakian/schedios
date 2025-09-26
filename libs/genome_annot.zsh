#!/usr/bin/env zsh

function annotate_genome_taxid() {
    # This script downloads a genome by Taxonomy ID from NCBI, softmasks it, and annotates genes using BRAKER3.
    # Dependencies: datasets, dataformat, jq, EDTA, BRAKER3

    # === Defaults ===
    threads=4
    prot_seq="${DATA}/OrthoDB11/viridiplantae.fa"

    # === Arguments ===
    if (( $# < 1 )); then
        echo "Usage: annotate_genome_taxid <taxid> [--prot_path <prot_path>] [--threads <threads>]" >&2
        exit 1
    fi
    taxid="$1"
    shift

    while (( $# > 0 )); do
        case "$1" in
            --threads)
                threads="$2"
                shift 2
                ;;
            --prot_path)
                prot_seq="$2"
                shift 2
                ;;
            *)
                echo "[Error] Unknown option: $1" >&2
                exit 1
                ;;
        esac
    done

    # === Obtain species name from Taxonomy ID ===
    json=$(datasets summary taxonomy taxon "$taxid")
    name=$(echo "$json" | jq -r '.reports[0].taxonomy.classification.species.name')
    if [[ -z "$name" || "$name" == "null" ]]; then
        echo "[ERROR] Failed to retrieve scientific name for taxid: $taxid" >&2
        exit 1
    fi
    sp_name="${name// /_}"

    # === Download genome data ===
    dir="${DATA}/${sp_name}"
    mkdir -p "$dir"
    datasets download genome taxon "$taxid" \
        --include genome \
        --reference \
        --filename "${dir}/${sp_name}.zip"

    # === Parse assembly data report ===
    dataformat tsv genome \
        --package "${dir}/${sp_name}.zip" \
        --fields organism-name,assminfo-name,accession \
        | tail -n +2 | IFS=$'\t' read -r org asm acc

    # === Process downloaded files ===
    unzip -o "${dir}/${sp_name}.zip" -d "$dir"
    datadir="${dir}/ncbi_dataset/data/${acc}"
    cp "${datadir}/${acc}_${asm}_genomic.fna" "${dir}/${sp_name}.dna.toplevel.fasta"

    # === Cleanup garbage ===
    rm -r "${dir}/ncbi_dataset"
    rm "${dir}/${sp_name}.zip"
    rm "${dir}/README.md"
    rm "${dir}/md5sum.txt"

    # === Softmask genome with RepeatMasker ===
    BuildDatabase -name "${dir}/database" \
        "${dir}/${sp_name}.dna.toplevel.fasta" \
        -pa "$threads"

    RepeatModeler -database "${dir}/database" \
        -pa "$(( threads / 4 ))"

    outdir=$(ls -td ${dir}/RM_* | head -n1)
    cat "${DATA}/Dfam/Dfam-RepeatMasker.lib" "${outdir}/consensi.fa.classified" > "${dir}/repeatmasker.lib"

    RepeatMasker -lib "${dir}/repeatmasker.lib" \
        -pa "$threads" \
        -xsmall \
        "${dir}/${sp_name}.dna.toplevel.fasta"

    # === Gene prediction with BRAKER3 ===
    workdir="${dir}/braker"
    braker3 \
        --genome="${dir}/${sp_name}.dna.toplevel.fasta.masked" \
        --prot_seq="$prot_seq" \
        --species="$sp_name" \
        --threads="$threads" \
        --workingdir="$workdir" \
        --gff3 \
        --AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH}"

    # === copy output files ===
    cp "${workdir}/braker.aa" "${dir}/${sp_name}.pep.all.fasta"
    cp "${workdir}/braker.codingseq" "${dir}/${sp_name}.cds.all.fasta"
    cp "${workdir}/braker.gff3" "${dir}/${sp_name}.genome.gff"
}