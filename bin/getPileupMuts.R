#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
nProcs <- as.numeric(args[3])
jsonFile <- as.character(args[4])
topPosFile <- as.character(args[5])
preMirFastaFile <- as.character(args[6])

source(file.path(scriptDir, "bin/functions.R"))

setupRlibs(R_libs)

library(dplyr)
library(purrr)
library(jsonlite)
library(readr)
library(tidyr)
library(Biostrings)
library(BiocParallel)

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

snow <- SnowParam(workers = nProcs, type = "SOCK")

topMirs <-
  topMirCutoff %>%
  group_by(flybase_id) %>%
  do(muts = pileupParallelMuts(groupedData = ., snow = snow)) %>%
  unnest(muts)

system.time(topMirs <- topMirCutoff %>% group_by(flybase_id) %>% do(muts = pileupParallelMuts(groupedData = ., snow = snow)) %>% unnest(muts))
