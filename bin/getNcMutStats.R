#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
nCores <- as.numeric(args[3])
jsonFile <- as.character(args[4])
posFile <- as.character(args[5])
fastaFile <- as.character(args[6])

source(file.path(scriptDir, 'bin/functions.R'))

# Packages loaded in setupRlibs

setupRlibs(R_libs)
pacman::p_load(Biostrings, Rsamtools, BiocParallel, jsonlite, tidyverse)
sessioninfo::session_info()

cfg <- getcfg(jsonFile)

mirPositions <- read_tsv(posFile)

fastaSeqs <- readDNAStringSet(fastaFile)
refTbl <- as_tibble(list('locus' = names(fastaSeqs),
	                     'full.seq' = paste(str_to_upper(fastaSeqs))))

agoSubstrMinLen <- 20
agoSubstrMaxLen <- 24

# Get tidy format of reference nucleotides from refTbl
maxSeqLen <- max(str_length(refTbl$full.seq))

# Create a tidy dataframe for one position per line, with its reference nucleotide
tidyRefNucs <-
  refTbl %>%
  separate(full.seq, paste('Pos', 1:maxSeqLen, sep = '_'), sep = '\\B') %>%
  gather(pos, refNuc, matches('Pos_')) %>% na.omit() %>%
  separate(pos, c('pos', 'idx'), sep = '_', convert = TRUE) %>%
  dplyr::select(-pos)

# Complete tibble with bamFiles, for effective do()
ncPosWFiles <-
  mirPositions %>%
  left_join(cfg$samples, by = c('timepoint' = 'id', 'time')) %>%
  mutate(bamFile = file.path(ncBam)) %>%
  dplyr::filter(!is.na(ncBam)) %>%
  left_join(refTbl) %>%
  dplyr::filter(average.ppm >= 5, read.type == "totalReads") %>%
  dplyr::rename(totalReads = reads)

# This is a bit hackey, as on our cluster, the nCores provided is half what R can use by hyperthreading
nProcs <- nCores * 2
mc.param <- MulticoreParam(workers = nProcs, type = 'FORK')

# Parallel pileup for all ncRNAs
ncsWmuts <-
  ncPosWFiles %>%
  group_by(bamFile) %>%
  do(muts = pileupMutsLenRestrict(groupedData = ., minLen = agoSubstrMinLen, maxLen = agoSubstrMaxLen, mc.param = mc.param)) %>%
  unnest(muts) %>%
  dplyr::select(-bamFile) %>%
  left_join(tidyRefNucs, by = c('locus', 'pos' = 'idx')) %>%
  left_join(ncPosWFiles %>% dplyr::select(-bamFile), by = c('locus', 'timepoint', 'time', 'start.pos' = 'pos'))

ncsWmuts %>% write_tsv('ncRNAs.wAllMuts.tsv')

# Determine mutation code and fraction of bases mutated
ncMutCodes <-
  ncsWmuts %>%
  spread(nucleotide, count) %>%
  replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
  gather(nucleotide, count, A:T) %>%
  mutate(mutCode = ifelse(refNuc != nucleotide,
                          paste(refNuc, nucleotide, sep = '>'),
                          refNuc)) %>%
  group_by(locus, timepoint, pos, start.pos) %>%
    mutate(depth = sum(count), mutFract = count / depth) %>%
    dplyr::filter(grepl('>', mutCode)) %>%
  ungroup() %>%
  dplyr::select(-refNuc, -nucleotide, -count, -align, -ncBam, -full.seq, -read.type)

# Switch to wide format for smoother excel copy/paste
ncMutsWide <-
  ncMutCodes %>%
  group_by(locus, time, start.pos, mutCode) %>%
    mutate(relMutCount = sprintf("%02d", rank(relPos)),
           relMut = paste(mutCode, relMutCount, sep = "_"),
           relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  dplyr::select(locus, start.pos, UCount, timepoint, time, average.ppm, depth, mutFract, relMut) %>%
  spread(relMut, mutFract) %>%
  arrange(locus, desc(average.ppm), time)

# Write out files
ncMutCodes %>% write_tsv('ncMutCodes.tsv')
ncMutsWide %>% write_tsv('ncMutsWide.tsv')
