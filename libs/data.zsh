#!/usr/bin/env zsh

function fetch_genome_by_genus() {
    # This script downloads genomes by genus from NCBI using the NCBI Datasets CLI.
    # Dependencies: datasets, biotp

    # === Usage ===
    if (( $# != 1 )); then
        echo "Usage: fetch_genome_by_genus <genus>" >&2
        exit 1
    fi

    # === Arguments ===
    genus="$1"

    # === Make directories ===
    mkdir -p ${DATA}/${genus}

    # === Download dataset ===
    datasets download genome taxon "$genus" \
        --include genome,protein,cds,gff3 \
        --reference \
        --annotated \
        --filename ${DATA}/${genus}/dataset.zip

    unzip -o ${DATA}/${genus}/dataset.zip -d ${DATA}/${genus}

    # == Parse assembly data report ==
    acc_li=()
    org_li=()
    asm_li=()

    tmpfile=$(mktemp)
    python3 -m biotp output_acc_org_asm \
        ${DATA}/${genus}/ncbi_dataset/data/assembly_data_report.jsonl > $tmpfile

    while IFS=" " read -r i acc org asm; do
        acc_li[i]="$acc" # acc: WGS accession
        org_li[i]="$org" # org: Scientific name
        asm_li[i]="$asm" # asm: Assembly
        echo "${i}, acc: ${acc_li[i]}, org: ${org_li[i]}, asm: ${asm_li[i]}"
    done < "$tmpfile"

    rm $tmpfile

    # === Copy files ===
    for ((i=1; i<=${#org_li[@]}; i++)); do
        acc="${acc_li[i]}"
        org="${org_li[i]}"
        asm="${asm_li[i]}"
        echo "[INFO] Processing organism: ${org} (Accession: ${acc}, Assembly: ${asm})"
        cp "${DATA}/${genus}/ncbi_dataset/data/${acc}/${acc}_${asm}_genomic.fna" \
            "${DATA}/${genus}/${org}.dna.toplevel.fasta"
        cp "${DATA}/${genus}/ncbi_dataset/data/${acc}/cds_from_genomic.fna" \
            "${DATA}/${genus}/${org}.cds.all.fasta"
        cp "${DATA}/${genus}/ncbi_dataset/data/${acc}/protein.faa" \
            "${DATA}/${genus}/${org}.pep.all.fasta"
        cp "${DATA}/${genus}/ncbi_dataset/data/${acc}/genomic.gff" \
            "${DATA}/${genus}/${org}.genome.gff"
    done
}
