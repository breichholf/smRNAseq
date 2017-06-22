# smRNAseq

Here we take a simple approach to map smallRNA-sequencing reads and map to miRNA hairpins downloaded form miRBase.

The pipeline is based on SciLifeLab/NGI-smRNAseq and borrows heavily from it. It allows passing of raw sequencing reads through the `--reads`, choosing a specific genome or version via `--genome` (currently only `dm6`), and passing a 0-3 mismatches using `--mismatches` [Default: 3].
