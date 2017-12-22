#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
jsonFile <- as.character(args[3])
preMirFastaFile <- as.character(args[4])

source(file.path(scriptDir, "bin/functions.R"))

setupRlibs(R_libs)

library(tidyverse)
library(Biostrings)
library(purrr)

sessionInfo()

cfg <- getcfg(jsonFile)

preMirFasta <- readDNAStringSet(preMirFastaFile)
preMirTbl <- as_tibble(list("flybase_id" = names(preMirFasta), "full.seq" = paste(preMirFasta)))

mirBodyLength <- 18

# `allcounts` records a count of all reads and T>C reads with T>C BQ>27 for all libraries.
# also includes lengths
allcounts <- cfg$samples %>% pmap(get.top.startpos, mirAnno = cfg$anno) %>% purrr::reduce(full_join)

# Converting to tidy format, and omitting length distribution.
gatheredCounts <-
  allcounts %>%
  select(-matches("LenDis")) %>%
  gather(type, reads, matches("Reads\\.")) %>%
  separate(type, c("read.type", "timepoint", "time"), sep = "\\.", convert = TRUE) %>%
  replace_na(list(reads = 0)) %>%
  group_by(flybase_id, pos, read.type, timepoint) %>%
    mutate(totalReads = sum(reads)) %>%
    select(-seqLen) %>% distinct() %>%
    filter(totalReads == max(totalReads)) %>%
  ungroup() %>% distinct()

# Get best mapping 5p starting position and add 'mature' and 'star' nomenclature
# If there are more than one position with equal reads, take the one closest to
# manually annotated 3p/5p start position (posdist)
topPositionCounts <-
  gatheredCounts %>%
  filter(read.type == "totalReads") %>%
  group_by(pos, flybase_id) %>%
    mutate(average.reads = mean(totalReads)) %>%
  group_by(arm.name) %>%
    top_n(1, average.reads) %>%
    mutate(posdist = ifelse(str_sub(arm.name, -2) == "3p",
                            `3p` - pos,
                            `5p` - pos)) %>%
  group_by(arm.name, average.reads) %>%
    filter(posdist == min(posdist)) %>%
    select(-posdist) %>%
  group_by(mir_name) %>%
    mutate(mir.type = ifelse(average.reads == max(average.reads), "mature", "star")) %>%
  ungroup()

# Adding in metadata such as U count and seed (pos 1-8)
topPosCntsWseed <-
  topPositionCounts %>%
  left_join(preMirTbl) %>%
  mutate(seed = str_sub(full.seq, pos, pos + 7),
         mirBody = str_sub(full.seq, pos, pos + mirBodyLength - 1),
         UCount = str_count(mirBody, "T")) %>%
  select(-mirBody, -full.seq, -reads) %>%
  distinct()

# Isolate TC reads
# topTcReads <-
#   gatheredCounts %>%
#   filter(read.type != "totalReads") %>%
#   left_join(topPosCntsWseed %>% select(flybase_id, pos, seed, UCount, timepoint, time, mir.type, average.reads)) %>%
#   filter(!is.na(mir.type))

# Convert `allcounts` to tidy format for length distributions (all reads and tc reads)
gatheredLenDis <-
  allCounts %>%
  select(-matches("Reads")) %>%
  gather(type, reads, matches("LenDis")) %>%
  separate(type, c("LD.type", "timepoint", "time"), sep = "\\.", convert = TRUE) %>%
  replace_na(list(reads = 0)) %>% distinct() %>%
  left_join(topPosCntsWseed %>% select(flybase_id, pos, seed, UCount, read.type, timepoint, time, mir.type, totalReads)) %>%
  filter(!is.na(mir.type))

# Isolate length distribution for TC reads only
# tcLenDis <-
#   gatheredLenDis %>%
#   filter(LD.type == "tcLenDis")

####
##  Write output files
####

# Separately extract mature and star steady state read counts (ppm) with metadata
matureWide <- convertToWide(topPosCntsWseed, "mature")
starWide <- convertToWide(topPosCntsWseed, "star")
# matureTcWide <- convertToWide(topTcReads, "mature")
# starTcWide <- convertToWide(topTcReads, "star")

allcounts %>% write_tsv('allCounts.tsv')
topPosCntsWseed %>% write_tsv('topPositionCounts.tsv')
gatheredLenDis %>% write_tsv('rawLenDis.tsv')
matureWide %>% write_tsv('miR.steadyState.PPM.tsv')
starWide %>% write_tsv('miRSTAR.steadyState.PPM.tsv')
# tcLenDis %>% write_tsv('tcLenDis.tsv')
# topTcReads %>% write_tsv('topTcReads.tsv')
# matureTcWide %>% write_tsv('matureTcPPM.tsv')
# starTcWide %>% write_tsv('starTcPPM.tsv')
