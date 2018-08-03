#!/usr/bin/env Rscript

possibleTypes <- c('rRNA', 'tRNA', 'snoRNA', 'snRNA', 'pre_miRNA')

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
lookupType <- as.character(args[3])
outName <- as.character(args[4])
downstreamNT <- if(length(args) == 5) as.numeric(args[5]) else 0

if (!(lookupType %in% possibleTypes)) {
    possTypeStr <- paste0("(", paste(possibleTypes, collapse = " "), ")")
    errorMessage <- paste("Lookup type", lookupType, "not amongst accepted types:", possTypeStr)
    stop(errorMessage)
}

source(file.path(scriptDir, "bin/functions.R"))

setupRlibs(R_libs)
pacman::p_load(biomaRt, seqinr, tidyverse)
sessioninfo::session_info()

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl")

attribList <- c('gene_biotype', 'ensembl_gene_id', 'external_gene_name', 'flybase_gene_id', 'flybase_transcript_id',
                'chromosome_name', 'start_position', 'end_position', 'strand')

martDf <- getBM(attributes = attribList, filters = 'biotype', values = lookupType, mart = ensembl)

seqDf <- biomaRt::getSequence(id = martDf$flybase_gene_id, type = "flybase_gene_id", downstream = downstreamNT,
                              seqType = "transcript_exon_intron", mart = ensembl)

martWseqs <-
  left_join(martDf, seqDf) %>%
  rename(seq = transcript_exon_intron) %>%
  mutate(longName = paste(flybase_gene_id, flybase_transcript_id, external_gene_name, gene_biotype, sep = "|"))

# We include the hackey end_position - start_position filter, as there is one snoRNA that misbehaves and is extracted twice
# One where the sequence length matches the difference between start and end, and another that doesn't
# Apart from that both entries are the same
# Using arrange(flybase_transcript_id) we ensure consistent naming across different runs of the script.
unique.seqs <-
  martWseqs %>%
  group_by(flybase_transcript_id) %>%
  filter(end_position - start_position + 1 + downstreamNT == str_length(seq)) %>%
  group_by(seq) %>%
  arrange(flybase_transcript_id) %>%
  filter(row_number() == min(row_number()))

write.fasta(as.list(unique.seqs$seq), unique.seqs$longName, outName)
