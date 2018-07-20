#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
nCores <- as.numeric(args[3])
jsonFile <- as.character(args[4])
posFile <- as.character(args[5])
preMirFastaFile <- as.character(args[6])

source(file.path(scriptDir, 'bin/functions.R'))

# Packages loaded in setupRlibs
#library(tidyverse)
#library(Biostrings)
#library(BiocParallel)

setupRlibs(R_libs)
pacman::p_load(Biostrings, Rsamtools, BiocParallel, jsonlite, tidyverse)
sessioninfo::session_info()

cfg <- getcfg(jsonFile)

mirPositions <- read_tsv(posFile)

preMirFasta <- readDNAStringSet(preMirFastaFile)
preMirTbl <- as_tibble(list('flybase_id' = names(preMirFasta), 'full.seq' = paste(preMirFasta)))

mirBodyLength <- 30

# Get tidy format of reference nucleotides from preMirTbl
maxHairpinLen <- max(str_length(preMirTbl$full.seq))

# Create a tidy dataframe for one position per line, with its reference nucleotide
tidyRefNucs <-
  preMirTbl %>%
  separate(full.seq, paste('Pos', 1:maxHairpinLen, sep = '_'), sep = '\\B') %>%
  gather(pos, refNuc, matches('Pos_')) %>% na.omit() %>%
  separate(pos, c('pos', 'idx'), sep = '_', convert = TRUE) %>%
  dplyr::select(-pos)

# Complete tibble with bamFiles, for effective do()
mirPosWFiles <-
  mirPositions %>%
  left_join(cfg$samples, by = c('timepoint' = 'id', 'time')) %>%
  mutate(bamFile = file.path(align)) %>%
  dplyr::filter(!is.na(align)) %>%
  left_join(preMirTbl) %>%
  dplyr::filter(average.ppm >= 5, read.type == "totalReads") %>%
  dplyr::rename(totalReads = reads)

# This is a bit hackey, as on our cluster, the nCores provided is half what R can use by hyperthreading
nProcs <- nCores * 2
mc.param <- MulticoreParam(workers = nProcs, type = 'FORK')

# Parallel pileup for all miRs
mirsWmuts <-
  mirPosWFiles %>%
  group_by(bamFile) %>%
  do(muts = pileupParallelMuts(groupedData = ., mc.param = mc.param, minLen = mirBodyLength)) %>%
  unnest(muts) %>%
  dplyr::select(-bamFile) %>%
  left_join(tidyRefNucs, by = c('flybase_id', 'pos' = 'idx')) %>%
  left_join(mirPosWFiles %>% dplyr::select(-bamFile), by = c('flybase_id', 'timepoint', 'time', 'mir.type', 'start.pos' = 'pos'))

mirsWmuts %>% write_tsv('miRs.wAllMuts.tsv')

# Determine mutation code and fraction of bases mutated
mirMutCodes <-
  mirsWmuts %>%
  spread(nucleotide, count) %>%
  replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
  gather(nucleotide, count, A:T) %>%
  mutate(mutCode = ifelse(refNuc != nucleotide,
                          paste(refNuc, nucleotide, sep = '>'),
                          refNuc)) %>%
  group_by(flybase_id, timepoint, pos, start.pos) %>%
    mutate(depth = sum(count), mutFract = count / depth) %>%
    dplyr::filter(grepl('>', mutCode)) %>%
  ungroup() %>%
  dplyr::select(-refNuc, -nucleotide, -count, -(`5p`:`3p`), -align, -ncBam, -full.seq, -mir_name, -read.type)

# Switch to wide format for smoother excel copy/paste
mirMutsWide <-
  mirMutCodes %>%
  group_by(flybase_id, time, start.pos, mutCode) %>%
    mutate(relMutCount = sprintf("%02d", rank(relPos)),
           relMut = paste(mutCode, relMutCount, sep = "_"),
           relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  dplyr::select(flybase_id, arm.name, mir.type, start.pos, seed, UCount, timepoint, time, average.ppm, depth, mutFract, relMut) %>%
  spread(relMut, mutFract) %>%
  arrange(mir.type, desc(average.ppm), time)

# Write out files
mirMutCodes %>% write_tsv('mirMutCodes.tsv')
mirMutsWide %>% write_tsv('mirMutsWide.tsv')
