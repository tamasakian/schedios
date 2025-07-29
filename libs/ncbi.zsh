#!/usr/bin/env zsh

function fetch_genome_taxid() {
    # This script downloads genomes by Taxonomy ID from NCBI using the NCBI Datasets CLI.
    # Dependencies: datasets, dataformat, jq

    # === Arguments ===
    if (( $# != 1 )); then
        echo "Usage: fetch_genome_taxid <taxid>" >&2
        exit 1
    fi
    taxid="$1"

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
        --include genome,protein,cds,gff3 \
        --reference \
        --annotated \
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
    cp "${datadir}/cds_from_genomic.fna"      "${dir}/${sp_name}.cds.all.fasta"
    cp "${datadir}/protein.faa"               "${dir}/${sp_name}.pep.all.fasta"
    cp "${datadir}/genomic.gff"               "${dir}/${sp_name}.genome.gff"

    # === Cleanup garbage ===
    rm -r "${dir}/ncbi_dataset"
    rm "${dir}/${sp_name}.zip"
    rm "${dir}/README.md"
    rm "${dir}/md5sum.txt"
}
