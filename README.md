# SCHEDIOS
SCHEDIOS: A Zsh-based sketching environment for computational craft in the life sciences.

## Usage
To use schedios, source the profile script in your shell:

```zsh
source ~/schedios/jobs/profile.zsh
```

### fetch_genome_by_genus
This function downloads genomes by genus from NCBI using the NCBI Datasets CLI.
It requires the [datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) and [biotp](https://github.com/tamasakian/biotp.git) commands to be available in your environment.
```zsh
fetch_genome_by_genus <genus>
```
- `<genus>`: The genus of the organism you want to fetch genomes for.
- The downloaded files will be renamed as follows:
  - `<genus>.<org>.dna.toplevel.fasta`: Genome sequence in FASTA format.
  - `<genus>.<org>.pep.all.fasta`: Protein sequences in FASTA format.
  - `<genus>.<org>.cds.all.fasta`: CDS sequences in FASTA format.
  - `<genus>.<org>.genome.gff`: Genome annotation in GFF format.

### forge_hmmscan_input
This function prepares input for HMMscan by creating a task directory and generating a combined FASTA file from multiple species.
It requires the [bithon](https://github.com/tamasakian/bithon.git) and [fasp](https://github.com/tamasakian/fasp.git) commands to be available in your environment.
```zsh
forge_hmmscan_input <sp_name1> <sp_name2> ...
```
- `<sp_name1>`, `<sp_name2>`, ...: Names of the species to include in the HMMscan input.
- The function creates a task directory with the current timestamp and generates a combined FASTA file named `all_species.fasta` in the `input` subdirectory of the task directory.
- The FASTA files for each species are prefixed with their species name to ensure unique identifiers.

