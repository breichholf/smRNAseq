#!/usr/bin/env Rscript

source("functions.R")

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
R_libs <- as.character(args[1])
jsonFile <- as.character(args[2])
gatherFile <- as.character(args[3])
preMirFastaFile <- as.character(args[4])

setupRlibs(R_libs)

library(dplyr)
library(purrr)
library(jsonlite)

cfg.info <- jsonlite::read_json(jsonFile)
file.home <- cfg.info$base
mir.anno <- read_tsv(cfg.info$mir.anno)

cfg.samples <- dplyr::bind_rows(cfg.info$samples)
cfg.samples <- mutate(cfg.samples, align = file.path(file.home, align))

gatheredCounts <- read_tsv(gatherFile)

preMirFasta <- readDNAStringSet(preMirFastaFile)
preMirTbl <- as_tibble(list("flybase_id" = names(preMirFasta), "full.seq" = paste(preMirFasta)))

topMirCounts <-
  gatheredCounts %>%
  left_join(cfg.samples, by = c('timepoint' = 'id')) %>%
  mutate(bamFile = file.path(file.home, align)) %>%
  group_by(arm.name) %>%
  top_n(1, average.reads) %>%
  ungroup() %>%
  left_join(preMirTbl)

topMirMuts <-
  topMirCounts %>%
  pmap_dfr(mutsFromPileup)

topMirMutCodes <-
  topMirMuts %>%
  spread(nucleotide, count) %>%
  replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
  gather(nucleotide, count, A:T) %>%
  mutate(mutCode = ifelse(refNuc != nucleotide,
                          paste(refNuc, nucleotide, sep = ">"),
                          refNuc)) %>%
  group_by(flybase_id, timepoint, pos) %>%
  mutate(depth = sum(count), mutFract = count / depth) %>%
  filter(grepl(">", mutCode)) %>%
  arrange(mutCode, relPos) %>% ungroup()

topMirMutCodes %>% write_tsv('mutstats.tsv')
