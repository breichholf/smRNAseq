#!/usr/bin/env Rscript

library(dplyr)
library(tibble)
library(purrr)
library(jsonlite)
library(Biostrings)

source("functions.R")

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
R_libs <- as.character(args[1])
jsonFile <- as.character(args[2])
topPosFile <- as.character(args[3])
preMirFile <- as.character(args[4])

setupRlibs(R_libs)

cfg.info <- jsonlite::read_json(jsonFile)
file.home <- cfg.info$base
mir.anno <- read_tsv(cfg.info$mir.anno)

cfg.samples <- dplyr::bind_rows(cfg.info$samples)
cfg.samples <- mutate(cfg.samples, align = file.path(file.home, align))

topPositions <- read_tsv(topPosFile)

preMirFasta <- readDNAStringSet(preMirFile)
preMirTbl <- as_tibble(list("flybase_id" = names(preMirFasta), "full.seq" = paste(preMirFasta)))

# !!!ATTENTION!!!
# This static filter does not take sRNA reads in to account!
# It should be cpm, but if we don't know sRNAreads, then this will reflect reads
topPosPastCutoff <- filter(topPositions, average.reads >= 50)

topPosWPath <-
  topPosPastCutoff %>%
  filter(!grepl("-2a-1|-2b-1", arm.name)) %>%
  left_join(cfg.samples, by = c("timepoint" = "id")) %>%
  mutate(bamFile = file.path(file.home, align)) %>%
  left_join(preMirTbl)

reducedSet <-
  topPosWPath %>%
  group_by(timepoint) %>%
  top_n(5, average.reads) %>% ungroup()

start.time <- Sys.time()
reducedSetWithMuts <- reducedSet %>% pmap_dfr(fetchPileup)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
