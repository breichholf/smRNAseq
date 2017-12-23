#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
allCountsFile <- as.character(args[3])
topPositionsFile <- as.character(args[4])

source(file.path(scriptDir, "bin/functions.R"))
setupRlibs(R_libs)

library(tidyverse)
sessionInfo()

mirBodyLength <- 18

# `allcounts` records a count of all reads and T>C reads with T>C BQ>27 for all libraries.
# also includes lengths
allCounts <- read_tsv(allCountsFile)
topPosCntsWseed <- read_tsv(topPositionsFile)

# Convert `allCounts` to tidy format for length distributions (all reads and tc reads)
gatheredLenDis <-
  allCounts %>%
  select(-matches("Reads")) %>%
  gather(type, reads, matches("LenDis")) %>%
  separate(type, c("LD.type", "timepoint", "time"), sep = "\\.", convert = TRUE) %>%
  replace_na(list(reads = 0)) %>% distinct() %>%
  left_join(topPosCntsWseed %>% select(flybase_id, pos, seed, UCount, read.type, timepoint, time, mir.type, reads, average.reads)) %>%
  filter(!is.na(mir.type))

####
##  Write output files
####

# Separately extract mature and star steady state read counts (ppm) with metadata
matureWide <- convertToWide(topPosCntsWseed, "mature")
starWide <- convertToWide(topPosCntsWseed, "star")

gatheredLenDis %>% write_tsv('rawLenDis.tsv')
matureWide %>% write_tsv('miR.steadyState.PPM.tsv')
starWide %>% write_tsv('miRSTAR.steadyState.PPM.tsv')
