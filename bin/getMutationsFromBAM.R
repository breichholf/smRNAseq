#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
jsonFile <- as.character(args[3])
topPosFile <- as.character(args[4])
preMirFastaFile <- as.character(args[5])

source(file.path(scriptDir, "bin/functions.R"))

setupRlibs(R_libs)

library(dplyr)
library(purrr)
library(jsonlite)
library(readr)
library(tidyr)
library(Biostrings)

sessionInfo()

cfg.info <- jsonlite::read_json(jsonFile)
file.home <- cfg.info$base
mir.anno <- read_tsv(cfg.info$mir.anno)

cfg.samples <- dplyr::bind_rows(cfg.info$samples)
cfg.samples <- mutate(cfg.samples, align = file.path(file.home, align))

topPositions <- read_tsv(topPosFile)

preMirFasta <- readDNAStringSet(preMirFastaFile)
preMirTbl <- as_tibble(list("flybase_id" = names(preMirFasta), "full.seq" = paste(preMirFasta)))

mirBodyLength <- 18

topMirCounts <-
  topPositions %>%
  left_join(cfg.samples, by = c('timepoint' = 'id')) %>%
  mutate(bamFile = file.path(align)) %>%
  left_join(preMirTbl)

topMirCutoff <-
  topMirCounts %>%
  dplyr::filter(average.reads >= 50)

topMirs <-
  topMirCutoff %>%
  pmap_dfr(mutsFromPileup, minLen = mirBodyLength) %>%
  dplyr::filter(relPos <= mirBodyLength)

# topMirMutsWarms <-
#   topMirMuts %>%
#   left_join(mir.anno, by = )

topMirMutCodes <-
  topMirs %>%
  spread(nucleotide, count) %>%
  replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
  gather(nucleotide, count, A:T) %>%
  mutate(mutCode = ifelse(refNuc != nucleotide,
                          paste(refNuc, nucleotide, sep = ">"),
                          refNuc)) %>%
  group_by(flybase_id, timepoint, pos) %>%
    mutate(depth = sum(count), mutFract = count / depth) %>%
    dplyr::filter(grepl(">", mutCode)) %>%
    arrange(flybase_id, timepoint, relPos, mutCode) %>%
  ungroup()

mirMutsWide <-
  topMirMutCodes %>%
  group_by(flybase_id, timepoint) %>%
    mutate(start.pos = min(pos)) %>%
  group_by(flybase_id, timepoint, mutCode) %>%
    mutate(relMutCount = rank(relPos)) %>%
    mutate(relMut = paste(mutCode, relMutCount, sep = "_"),
           relMut = str_replace(relMut, ">", "")) %>% ungroup() %>%
    dplyr::select(flybase_id, timepoint, start.pos, depth, mutFract, relMut) %>%
  spread(relMut, mutFract)

topMirMutCodes %>% write_tsv('mutStats.tsv')
mirMutsWide %>% write_tsv('allMirMuts.tsv')
