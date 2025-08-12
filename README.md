# SCHEDIOS
SCHEDIOS: A Zsh-based sketching environment for computational craft in the life sciences.

## Usage
To use schedios, source the profile script in your shell:

```zsh
source ~/schedios/jobs/profile.zsh
```

### fetch_genome_taxid
This function fetches the NCBI Taxonomy ID (TaxID) for a given genome name.
It requires the `datasets` and `dataformat` command to be available in your environment.
```zsh
fetch_genome_taxid taxid
```
- `taxid`: The NCBI Taxonomy ID of the genome you want to fetch.
- The downloaded files will be renamed as follows:
  - `dna.genome.fasta`: The genome sequence in FASTA format.
  - `genome.gff`: The genome annotation in GFF format.
  - `cds.all.fasta`: The coding sequences in FASTA format.
  - `pep.all.fasta`: The protein sequences in FASTA format.

### run_orthofinder
This function runs OrthoFinder on a set of genomes.
It requires the `bython` command to be available in your environment.

```zsh
input4orthofinder sp_name1 sp_name2 sp_name3 ...
```
- `sp_name1`, `sp_name2`, `sp_name3`, ...: The names of the species to include in the OrthoFinder analysis.
- The function will create a task directory named `orthofinder` and download the genomes for the specified species.
- The genomes will be stored in the `orthofinder/input` directory.

```zsh
run_orthofinder task_name
```
- `task_name`: The name of the task directory where the OrthoFinder results will be stored.
- The function will run OrthoFinder and store the results in the specified task directory.