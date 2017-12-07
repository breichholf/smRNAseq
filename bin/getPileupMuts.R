#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
nProcs <- as.numeric(args[3])
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
  left_join(preMirTbl)

topMirCutoff <-
  topMirCounts %>%
  dplyr::filter(average.reads >= 50)

snow <- SnowParam(workers = nProcs, type = 'SOCK')

topMirs <-
  topMirCutoff %>%
  group_by(bamFile) %>%
  do(muts = pileupParallelMuts(groupedData = ., snow = snow)) %>%
  unnest(muts) %>%
  dplyr::select(-bamFile) %>%
  left_join(tidyRefNucs, by = c('flybase_id', 'pos' = 'idx'))
  # left_join(topMirCutoff, by = c('flybase_id', 'timepoint', 'mir.type', 'start.pos' = 'pos'))

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
    arrange(flybase_id, time, relPos, mutCode) %>%
  ungroup()

mirMutsWide <-
  topMirMutCodes %>%
  group_by(flybase_id, timepoint, mutCode) %>%
    mutate(relMutCount = sprintf("%02d", rank(relPos)),
           relMut = paste(mutCode, relMutCount, sep = "_"),
           relMut = str_replace(relMut, ">", "")) %>% ungroup() %>%
    dplyr::select(flybase_id, timepoint, time, start.pos, depth, mutFract, relMut) %>%
  spread(relMut, mutFract) %>%
  arrange(flybase_id, start.pos, time)

topMirMutCodes %>% write_tsv('mutStats.tsv')
mirMutsWide %>% write_tsv('allMirMuts.tsv')
