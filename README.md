# SCHEDIOS
SCHEDIOS: A Zsh-based sketching environment for computational craft in the life sciences.

# Usage
To use schedios, source the profile script in your shell:

```zsh
source ~/schedios/jobs/profile.zsh
```

## fetch_genome_by_genus
This function downloads genomes by genus from NCBI using the NCBI Datasets CLI.
It requires the `datasets` and `biotp` commands to be available in your environment.
```zsh
fetch_genome_by_genus <genus>
```
- `<genus>`: The genus of the organism you want to fetch genomes for.
- The downloaded files will be renamed as follows:
  - `<genus>.<org>.dna.toplevel.fasta`: Genome sequence in FASTA format.
  - `<genus>.<org>.pep.all.fasta`: Protein sequences in FASTA format.
  - `<genus>.<org>.cds.all.fasta`: CDS sequences in FASTA format.
  - `<genus>.<org>.genome.gff`: Genome annotation in GFF format.

