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

cfg.info <- jsonlite::read_json(jsonFile)
file.home <- cfg.info$base
mir.anno <- read_tsv(cfg.info$mir.anno)

cfg.samples <- dplyr::bind_rows(cfg.info$samples)
cfg.samples <- mutate(cfg.samples, align = file.path(file.home, align))

allcounts <- cfg.samples %>% pmap(calc.maxpos, mirAnno = mir.anno) %>% purrr::reduce(full_join)

gatheredCounts <-
  allcounts %>%
  gather(type, reads, matches("Reads\\.")) %>%
  separate(type, c("read.type", "timepoint"), convert = TRUE) %>%
  replace_na(list(reads = 0)) %>%
  group_by(pos, flybase_id, read.type) %>%
  mutate(average.reads = mean(reads))

topPositionCounts <-
  gatheredCounts %>%
  group_by(arm.name) %>% top_n(1, average.reads) %>%
  group_by(mir_name) %>%
  mutate(mir.type = ifelse(average.reads == min(average.reads), "star", "mature")) %>%
  ungroup()

matureWide <- convertToWide(topPositionCounts, "mature")

starWide <- convertToWide(topPositionCounts, "star")

allcounts %>% write_tsv('allCounts.tsv')
gatheredCounts %>% write_tsv('gatheredCounts.tsv')
topPositionCounts %>% write_tsv('topPositions.tsv')
matureWide %>% write_tsv('matureMirs.tsv')
starWide %>% write_tsv('starMirs.tsv')
