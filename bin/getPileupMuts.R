#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
nCores <- as.numeric(args[3])
jsonFile <- as.character(args[4])
topPosFile <- as.character(args[5])
preMirFastaFile <- as.character(args[6])

source(file.path(scriptDir, 'bin/functions.R'))

setupRlibs(R_libs)

library(tidyverse)
library(Biostrings)
library(BiocParallel)

sessionInfo()

cfg <- getcfg(jsonFile)

topPositions <- read_tsv(topPosFile)

preMirFasta <- readDNAStringSet(preMirFastaFile)
preMirTbl <- as_tibble(list('flybase_id' = names(preMirFasta), 'full.seq' = paste(preMirFasta)))

# Get tidy format of reference nucleotides from preMirTbl
maxHairpinLen <- max(str_length(preMirTbl$full.seq))

tidyRefNucs <-
  preMirTbl %>%
  separate(full.seq, paste('Pos', 1:maxHairpinLen, sep = '_'), sep = '\\B') %>%
  gather(pos, refNuc, matches('Pos_')) %>% na.omit() %>%
  separate(pos, c('pos', 'idx'), sep = '_', convert = TRUE) %>%
  dplyr::select(-pos)

topMirCounts <-
  topPositions %>%
  left_join(cfg$samples, by = c('timepoint' = 'id', 'time')) %>%
  mutate(bamFile = file.path(align)) %>%
  filter(!is.na(align)) %>%
  left_join(preMirTbl)

topMirCutoff <-
  topMirCounts %>%
  dplyr::filter(average.reads >= 50)

nProcs <- nCores * 2
mc.param <- MulticoreParam(workers = nProcs, type = 'FORK')

topMirs <-
  topMirCutoff %>%
  group_by(bamFile) %>%
  do(muts = pileupParallelMuts(groupedData = ., mc.param = mc.param)) %>%
  unnest(muts) %>%
  dplyr::select(-bamFile) %>%
  left_join(tidyRefNucs, by = c('flybase_id', 'pos' = 'idx')) %>%
  left_join(topMirCutoff %>% select(-bamFile), by = c('flybase_id', 'timepoint', 'time', 'mir.type', 'start.pos' = 'pos'))

topMirMutCodes <-
  topMirs %>%
  spread(nucleotide, count) %>%
  replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
  gather(nucleotide, count, A:T) %>%
  mutate(mutCode = ifelse(refNuc != nucleotide,
                          paste(refNuc, nucleotide, sep = '>'),
                          refNuc)) %>%
  group_by(flybase_id, timepoint, pos) %>%
    mutate(depth = sum(count), mutFract = count / depth) %>%
    dplyr::filter(grepl('>', mutCode)) %>%
  ungroup() %>%
  select(-refNuc, -nucleotide, -count, -`5p`, -`3p`, -align, -full.seq, -mir_name, -read.type)

mirMutsWide <-
  topMirMutCodes %>%
  group_by(flybase_id, time, start.pos, mutCode) %>%
    mutate(relMutCount = sprintf("%02d", rank(relPos)),
           relMut = paste(mutCode, relMutCount, sep = "_"),
           relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  dplyr::select(flybase_id, arm.name, mir.type, start.pos, seed, UCount, timepoint, time, average.reads, depth, mutFract, relMut) %>%
  spread(relMut, mutFract) %>%
  arrange(mir.type, desc(average.reads), time)

topMirMutCodes %>% write_tsv('mutStats.tsv')
mirMutsWide %>% write_tsv('allMirMuts.tsv')
