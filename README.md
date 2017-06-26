# smRNAseq

Here we take a simple approach to map smallRNA-sequencing reads and map to miRNA hairpins downloaded form miRBase.

The pipeline is based on SciLifeLab/NGI-smRNAseq and borrows heavily from it. It allows passing of raw sequencing reads through the `--reads`, choosing a specific genome or version via `--genome` (currently only `dm6`), and passing a 0-3 mismatches using `--mismatches` [Default: 3].

# Planned features

In no particular order:

[x] Write nextflow pipeline for `bowtie -> arrayTagTCreads.awk -> write json`
[x] JSON for file parsing should give all information
[ ] Get mapping statistics with QuasR(?)
[ ] Get U-list for all positions from 0 - END-21 (python -- script ready)
[ ] Parse JSON and create count & TC counts
[ ] Get U>C conversion (per position) from pileup
[ ] Get all mutations from pileups
[ ] Tally all mutations in boxplot
[ ] Combining R analysis and nextflow: json format for files?
[ ] Extract smallRNA mappers from custom bash-AnnotationPipeline?
