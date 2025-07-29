#!/usr/bin/env zsh

function fetch_genome_taxid() {
    # This script downloads genomes by Taxonomy ID from NCBI using the NCBI Datasets CLI.
    # Dependencies: datasets, dataformat
    if (( $# != 1 )); then
        echo "Usage: fetch_genome_taxid <Taxonomy ID>" >&2
        exit 1
    fi

    # === Arguments ===
    taxid="$1"

    # === Obtain species and genus names from Taxonomy ID ===
    json=$(datasets summary taxonomy taxon "$taxid")
    name=$(echo "$json" | jq -r '.reports[0].taxonomy.current_scientific_name.name')
    if [[ -z "$name" || "$name" == "null" ]]; then
        echo "[ERROR] Failed to retrieve scientific name for taxid: $taxid" >&2
        exit 1
    fi
    gn_name="${name%% *}"
    sp_name="${name// /_}"

    echo "[INFO] Taxonomy ID: $taxid"
    echo "[INFO] Scientific name: $name"
    echo "[INFO] Genus name: $gn_name"
    echo "[INFO] Species name: $sp_name"

    # === Make directories ===
    dir="${DATA}/${gn_name}"
    mkdir -p "$dir"

    # === Download ===
    datasets download genome taxon "$taxid" \
        --include genome,protein,cds,gff3 \
        --reference \
        --annotated \
        --filename "${dir}/${sp_name}.zip"

    unzip -o "${dir}/${sp_name}.zip" -d "$dir"

    # === Parse assembly data report ===
    dataformat tsv genome \
        --package "${dir}/${sp_name}.zip" \
        --fields organism-name,assminfo-name,accession \
        | tail -n +2 | IFS=$'\t' read -r org asm acc

    echo "[INFO] Assembly: $asm"
    echo "[INFO] Accession: $acc"

    datadir="${dir}/ncbi_dataset/data/${acc}"
    mv "${datadir}/${acc}_${asm}_genomic.fna" "${dir}/${sp_name}.dna.toplevel.fasta"
    mv "${datadir}/cds_from_genomic.fna"      "${dir}/${sp_name}.cds.all.fasta"
    mv "${datadir}/protein.faa"               "${dir}/${sp_name}.pep.all.fasta"
    mv "${datadir}/genomic.gff"               "${dir}/${sp_name}.genome.gff"

    # === Cleanup ===
    rm -rf "${dir}/ncbi_dataset"
    rm "${dir}/${sp_name}.zip"
}
