#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
jsonFile <- as.character(args[3])
gatherFile <- as.character(args[4])
preMirFastaFile <- as.character(args[5])

source(file.path(scriptDir, "bin/functions.R"))

setupRlibs(R_libs)

library(dplyr)
library(purrr)
library(jsonlite)
library(readr)

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

# topMirMutsWarms <-
#   topMirMuts %>%
#   left_join(mir.anno, by = )

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
  arrange(flybase_id, timepoint, relPos, mutCode) %>% ungroup()

mirMutsWide <-
  topMirMutCodes %>%
  group_by(flybase_id, timepoint) %>%
  mutate(start.pos = min(pos)) %>%
  group_by(flybase_id, timepoint, mutCode) %>%
  mutate(relMutCount = rank(relPos)) %>%
  mutate(relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", "")) %>% ungroup() %>%
  select(flybase_id, timepoint, start.pos, depth, mutFract, relMut) %>%
  spread(relMut, mutFract)

topMirMutCodes %>% write_tsv('mutstats.tsv')
mirMutsWide %>% write_tsv('allMirMuts.tsv')
