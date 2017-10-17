#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
jsonFile <- as.character(args[3])

source(file.path(scriptDir, "bin/functions.R"))

setupRlibs(R_libs)

library(dplyr)
library(purrr)
library(jsonlite)
library(tidyr)
library(readr)

sessionInfo()

cfg.info <- jsonlite::read_json(jsonFile)
file.home <- cfg.info$base
mir.anno <- read_tsv(cfg.info$mir.anno)

cfg.samples <- dplyr::bind_rows(cfg.info$samples)
cfg.samples <- mutate(cfg.samples, align = file.path(file.home, align))

allcounts <- cfg.samples %>% pmap(get.top.startpos, mirAnno = mir.anno) %>% purrr::reduce(full_join)

gatheredCounts <-
  allcounts %>%
  select(-matches("LenDis"), -seqLen) %>%
  gather(type, reads, matches("Reads\\.")) %>%
  separate(type, c("read.type", "timepoint"), convert = TRUE) %>%
  replace_na(list(reads = 0)) %>%
  group_by(flybase_id, pos, read.type, timepoint) %>%
    dplyr::filter(reads == max(reads)) %>%
  ungroup() %>% distinct()

topPositionCounts <-
  gatheredCounts %>%
  dplyr::filter(read.type == "totalReads") %>%
  group_by(pos, flybase_id) %>%
    mutate(average.reads = mean(reads)) %>%
  group_by(arm.name) %>%
    top_n(1, average.reads) %>%
  group_by(mir_name) %>%
    mutate(mir.type = ifelse(average.reads == max(average.reads), "mature", "star")) %>%
  ungroup()

topTcReads <-
  gatheredCounts %>%
  dplyr::filter(read.type != "totalReads") %>%
  left_join(topPositionCounts %>% select(flybase_id, pos, mir.type)) %>%
  dplyr::filter(!is.na(mir.type))

gatheredLenDis <-
  allcounts %>%
  select(-matches("Reads")) %>%
  gather(type, reads, matches("LenDis")) %>%
  separate(type, c("LD.type", "timepoint"), convert = TRUE) %>%
  replace_na(list(reads = 0)) %>% distinct() %>%
  left_join(topPositionCounts %>% select(flybase_id, pos, mir.type)) %>%
  dplyr::filter(!is.na(mir.type))

totalLenDis <-
  gatheredLenDis %>%
  dplyr::filter(LD.type == "totalLenDis")

tcLenDis <-
  gatheredLenDis %>%
  dplyr::filter(LD.type == "tcLenDis")

# matureWide <- convertToWide(topPositionCounts, "mature")

# starWide <- convertToWide(topPositionCounts, "star")

allcounts %>% write_tsv('allCounts.tsv')
gatheredCounts %>% write_tsv('gatheredCounts.tsv')
totalLenDis %>% write_tsv('totalLenDis.tsv')
tcLenDis %>% write_tsv('tcLenDis.tsv')
topPositionCounts %>% write_tsv('topPositionCounts.tsv')
topTcReads %>% write_tsv('topTcReads.tsv')
# matureWide %>% write_tsv('matureMirs.tsv')
# starWide %>% write_tsv('starMirs.tsv')
