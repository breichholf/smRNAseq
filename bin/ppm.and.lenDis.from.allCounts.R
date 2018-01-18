#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
allCountsFile <- as.character(args[3])
preMirFastaFile <- as.character(args[4])

source(file.path(scriptDir, "bin/functions.R"))
setupRlibs(R_libs)

library(tidyverse)
library(Biostrings)
sessionInfo()

preMirFasta <- readDNAStringSet(preMirFastaFile)
preMirTbl <- as_tibble(list("flybase_id" = names(preMirFasta), "full.seq" = paste(preMirFasta)))

mirBodyLength <- 18

# `allcounts` records a count (in ppm normalised to smallRNAs) of:
# (1) all reads and (2) T>C reads with T>C BQ>27 for all libraries.
# It lists all sequences starting at `pos`, with length of `seqLen`,
# mapping to `flybase_id`, alias `mir_name`.
# mir arms are named `-3p` and `-5p` according to the manual annotation from mirbase
# provided in $mirArmAnno (see main.nf)
allCounts <- read_tsv(allCountsFile)

# Columns with `Reads` in them already contain summarised read counts
cleanAllCounts <- allCounts %>% select(-`5p`, -`3p`, -matches("Reads"))

gatheredCounts <-
  cleanAllCounts %>%
  gather(type, ppm, -pos, -seqLen, -flybase_id, -mir_name, -arm.name) %>%
  replace_na(list(ppm = 0)) %>%
  separate(type, c("count.type", "timepoint", "time"), sep = "\\.", convert = TRUE)

# `read.sum.wMetaInfo` is of tidy format and contains cleaned `read.sum`
# Added metadata consists of:
#   (1) `seed`, starting at `pos`
#   (2) `ucount`, counted over `mirBodyLength` (18nt) starting at `pos`
read.sum.wMetaInfo <-
  gatheredCounts %>%
  group_by(pos, flybase_id, count.type, timepoint) %>%
    mutate(read.sum = sum(ppm)) %>%
  ungroup() %>%
  select(-seqLen, -ppm) %>% distinct() %>%
  left_join(preMirTbl) %>%
  mutate(mirBody = str_sub(full.seq, pos, pos + mirBodyLength - 1),
         seed = str_sub(mirBody, 1, 8),
         ucount = str_count(mirBody, "T")) %>%
  select(-full.seq, -mirBody)

# Save gatheredCounts and read.sum.wMetaInfo as raw output first
# Further processing will be done in a separate file
