# smRNAseq - reporting mutations across the most abundant miRNA species

Inspired by SciLifeLab/NGI-smRNAseq.

Here we take a simple approach to map smallRNA-sequencing reads and map to miRNA hairpins downloaded form miRBase.

The following command line parameters can be used to customise running of the pipeline:
- `--reads` ... provide the read files to use as input. Wildcard possible, but neads to be provided in single quotation marks `[Required]`
- `--genome` ... specify the genome (or version) to use. Currently only compatible with `dm6` `[Required]`
- `--mismatches` ... define the number of mismatches to use (Possible values: 0-3, Default: 3)
- `--adapter` ... Adapter sequence to use for clipping (Default: Illumina TrueSeq 3' adapter + index primer: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG`)
- `--norand` ... Use this to **prevent** trimming of `NNNN` from both 5' and 3' end of reads
- `--annoreads` ... Sample file, providing columns `idx`, `sRNAreads`, `miRNAreads`, `time` for each sample `[Required]`

## Input file format

The pipeline assumes a certain file format: all input files are labelled `IDX_filename.fq` (`.gz` files are handelled by the pipeline transparently). This same IDX needs to separated from the remaining file name by an `_` (underscore) and is expected in the `idx` column of the sample file passed through `--annoreads`.

## Sample file format (passed through `--annoreads`)

The sample file requires the following columns (in no specific order): `idx`, `sRNAreads`, `miRNAreads`, `time`. To identify each sample uniquely, `idx` needs to have a unique value, and be matched with the `IDX` section in the input file names. The values for `sRNAreads` and `miRNAreads` can be pre-calculated by external tools. The `time` column is helpful for downstream analysis, if analysing time course experiments, but currently not evaluated in the scope of this pipeline.

## Read normalisation (within pipeline)

### `sRNAreads`

Currently for fly data, we map our smallRNA reads to the genome, perform hirarchical assignment of their annotation class (e.g. rRNA, tRNA, mitochondiral, snoRNA, snRNA, miRNA, piRNA, exonic, intronic, TEs) and then use the sum of all smallRNAs mapping reads as normalisation factor.

### `miRNAreads`

Currently not used, but provided for optional normalisation in downstream exploratory analyses.

## Example file formatting

Given three samples, that were in a time course and taken after `0`, `1` and `2` hours of treatment, the sequencing data could be labelled as follows:

Input files: `1234_sample1_0h.fq`, `1235_sample2_1h.fq`, `1236_sample3_2h.fq`

After mapping and categorising all small RNA and miRNA mapping reads, we could have the following summary table:

|      filename      | sRNAreads | miRNAreads |
|--------------------|-----------|------------|
| 1234_sample1_0h.fq |  10239489 |    8172363 |
| 1235_sample2_1h.fq |  14864299 |   10236982 |
| 1236_sample3_2h.fq |  11928376 |    9283726 |
|--------------------|-----------|------------|

The manually created sample file, passed after the `--annoreads` parameter should therefore be formatted as follows:

| idx  | sRNAreads | miRNAreads | time |
|------|-----------|------------|------|
| 1234 |  10239489 |    8172363 |    0 |
| 1235 |  14864299 |   10236982 |    1 |
| 1236 |  11928376 |    9283726 |    2 |

# Planned features

- Compatibility
    + [ ] Add other genomes (specifically: mouse & human)
- Preprocessing
    + [ ] Nextflow sequential annotation pipeline
        1. ribosomal DNA
        2. Mitochondrial genome
        3. genome unique, then multimappers
        4. pre-miRs unique, multimappers
    + [ ] Quantify from above -- extract smallRNA reads for normalisation
    + [ ] Extract smallRNA mappers from custom bash-AnnotationPipeline?
- Analysis/Reporting
    + [ ] Tally all mutations in boxplot
    + [ ] Combining R analysis and nextflow
- Documentation
    + [ ] Installation of the pipeline
    + [ ] Adding custom genomes
- Usability
    + [ ] Unit tests?
