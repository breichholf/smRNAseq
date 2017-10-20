#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
allCntFile <- as.character(args[1])
topPosFile <- as.character(args[2])
outDir <- as.character(args[3])

library(tidyverse)

topPosMirs <- read_tsv(topPosFile)

topPosMirs %>% unite(lib, read.type, timepoint, sep = ".") %>%
  spread(lib, reads) %>%
  write_tsv(file.path(outDir, 'spreadTopPos.tsv'))

allCounts <- read_tsv(allCntFile)

allLD <-
  allCounts %>% select(-matches('Reads\\.')) %>%
  gather(type, reads, matches("LenDis")) %>%
  separate(type, c("LD.type", "timepoint"), convert = TRUE) %>%
  left_join(topPosMirs %>% select(flybase_id, pos, timepoint, mir.type, average.reads) %>% distinct()) %>%
  filter(!is.na(mir.type))

bgLD <-
  allLD %>% filter(timepoint == min(timepoint)) %>%
  select(pos, seqLen, flybase_id, LD.type, reads) %>%
  replace_na(list(reads = 0))

allReadsBGMinus <-
  allLD %>% left_join(bgLD %>% rename(bg.reads = reads)) %>%
  mutate(bg.subtract = ifelse(reads - bg.reads > 0, reads - bg.reads, 0)) %>%
  group_by(pos, flybase_id, LD.type, timepoint) %>%
  mutate(read.sum = sum(bg.subtract, na.rm = TRUE)) %>%
  select(-seqLen, -reads, -bg.reads, -bg.subtract) %>%
  distinct() %>% ungroup()

maxReads <-
  allReadsBGMinus %>%
  filter(timepoint == max(timepoint), mir.type == "mature") %>%
  select(pos, flybase_id, LD.type, read.sum) %>%
  rename(max.reads = read.sum)

allReadsBGMinus %>%
  filter(LD.type != "totalLenDis") %>%
  unite(lendis, LD.type, timepoint, sep = ".") %>%
  spread(lendis, read.sum) %>%
  write_tsv(file.path(outDir, 'bgSubtractedReads.tsv'))

allReadsBGMinus %>%
  filter(LD.type != "totalLenDis") %>%
  left_join(maxReads) %>%
  mutate(read.norm = read.sum / max.reads) %>% select(-read.sum, max.reads) %>%
  unite(lendis, LD.type, timepoint, sep = ".") %>%
  spread(lendis, read.norm) %>%
  write_tsv(file.path(outDir, 'bgSubtractedNormReads.tsv'))
