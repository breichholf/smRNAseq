#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
jsonFile <- as.character(args[3])
fastaFile <- as.character(args[4])
minLen <- as.numeric(args[5])
maxLen <- as.numeric(args[6])

source(file.path(scriptDir, "bin/functions.R"))

# Packages loaded in setupRlibs
# library(tidyverse)
# library(Biostrings)

setupRlibs(R_libs)
pacman::p_load(Biostrings, Rsamtools, BiocParallel, jsonlite, tidyverse)
sessioninfo::session_info()

cfg <- getcfg(jsonFile)

fastaSeqs <- readDNAStringSet(fastaFile)
refTbl <- as_tibble(list("locus" = names(fastaSeqs),
                         "full.seq" = paste(str_to_upper(fastaSeqs))))

seqBodyLength <- 18

# How many IDs do we have? Should correspond with unique entries in cfg$samples
ids <- cfg$samples %>% select(id) %>% distinct()
idCount <- dim(ids)[1]

# `allcounts` records a count of all reads and T>C reads with T>C BQ>27 for all libraries.
# General columns (pos = start position):
#  - pos, seqLen, locus
#
# Columns (per ID) -- all normalised to sRNAreads from cfg$samples:
#  - tcLenDis.<libID>    tcReads.<libID>    totalLenDis.<libID>    totalReads.<libID>
allCounts <- cfg$samples %>% pmap(getNcRNACounts, minLen = minLen, maxLen = maxLen) %>% purrr::reduce(full_join)

allCounts %>% write_tsv('ncRNACounts.tsv')

# Long format of tcReads and totalReads - only 1 entry per starting position and miR name
# In the wide format, there might be some lengths that are present in one timepoint but not the other.
# When converting to long, these would be represented as NA in 'totalReads'. To get around that,
# we first convert NA to 0, select all unique entries and then sum up all reads.
gatheredCounts <-
  allCounts %>%
  select(-matches('LenDis')) %>%
  gather(type, reads, matches("Reads\\.")) %>%
  replace_na(list(reads = 0)) %>%
  separate(type, c("read.type", "timepoint", "time"), sep = "\\.", convert = TRUE) %>%
  group_by(pos, locus, read.type, timepoint) %>%
    select(-seqLen) %>% distinct() %>%
    mutate(readSum = sum(reads)) %>% select(-reads) %>% distinct() %>%
    dplyr::rename(reads = readSum) %>%
  ungroup()

# Final filtering: we're only interested in miRs that have > 0 at any timepoint.
nonZeroPos <-
  gatheredCounts %>%
  filter(read.type == "totalReads", reads != 0) %>%
  group_by(pos, locus) %>% filter(n() == idCount) %>% ungroup() %>%
  select(pos, locus, reads) %>% distinct() %>%
  group_by(pos, locus) %>% mutate(average.ppm = mean(reads)) %>%
  distinct() %>%
  left_join(refTbl) %>%
  mutate(seqBody = str_sub(full.seq, pos, pos + seqBodyLength - 1),
         UCount = str_count(seqBody, "T")) %>%
  select(-seqBody, -full.seq)

# Last but not least we add metadata (seed & UCount) for each starting position
gatheredCounts.gtZero <-
  nonZeroPos %>%
  left_join(gatheredCounts) %>%
  distinct()

# Filter length distributions
gatheredLenDis <-
  allCounts %>%
  select(-matches('Reads')) %>%
  gather(type, reads, matches("LenDis")) %>%
  replace_na(list(reads = 0)) %>% distinct() %>%
  separate(type, c("LD.type", "timepoint", "time"), sep = "\\.", convert = TRUE)

gatheredLenDis.gtZero <-
  nonZeroPos %>%
  left_join(gatheredLenDis) %>%
  distinct()

# Save all files -- there is no label yet which is mature and which is star!
gatheredCounts.gtZero %>% write_tsv('filteredNcPositions.tsv')
gatheredLenDis.gtZero %>% write_tsv('filteredNcLenDis.tsv')
