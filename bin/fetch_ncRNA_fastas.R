#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])

source(file.path(scriptDir, "bin/functions.R"))

setupRlibs(R_libs)
pacman::p_load(biomaRt, seqinr, tidyverse)
sessioninfo::session_info()

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl")

lookupTypes <- c('tRNA', 'rRNA', 'snoRNA', 'snRNA')

attribList <- c('gene_biotype', 'ensembl_gene_id', 'external_gene_name', 'flybase_gene_id', 'flybase_transcript_id',
                'chromosome_name', 'start_position', 'end_position', 'strand')

martDf <- getBM(attributes = attribList, filters = 'biotype', values = lookupTypes, mart = ensembl)

seqDf <- biomaRt::getSequence(id = martDf$flybase_gene_id, downstream = 20,
    type = "flybase_gene_id", seqType = "transcript_exon_intron", mart = ensembl)

martWseqs <-
  left_join(martDf, seqDf) %>%
  rename(seq = transcript_exon_intron) %>%
  mutate(longName = paste(flybase_gene_id, flybase_transcript_id, external_gene_name, gene_biotype, sep = "|"))

rRNA.unique <- martWseqs %>% filter(gene_biotype == "rRNA") %>% group_by(seq) %>% top_n(1)
tRNA.unique <- martWseqs %>% filter(gene_biotype == "tRNA") %>% group_by(seq) %>% top_n(1)
snRNA.unique <- martWseqs %>% filter(gene_biotype == "snRNA") %>% group_by(seq) %>% top_n(1)
snoRNA.unique <- martWseqs %>% filter(gene_biotype == "snoRNA") %>% group_by(seq) %>% top_n(1) %>%
  group_by(flybase_transcript_id) %>% filter(end_position - start_position + 21 == str_length(seq))

write.fasta(as.list(rRNA.unique$seq), rRNA.unique$longName, "ribosomes.fa")
write.fasta(as.list(tRNA.unique$seq), tRNA.unique$longName, "tRNA.fa")
write.fasta(as.list(snRNA.unique$seq), snRNA.unique$longName, "snRNA.fa")
write.fasta(as.list(snoRNA.unique$seq), snoRNA.unique$longName, "snoRNA.fa")
