#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
rawTcLenDisFile <- as.character(args[2])
topPosFile <- as.character(args[3])
normLDtime <- as.numeric(args[4])
maxTime <- as.numeric(args[5])
outDir <- as.character(args[6])

source(file.path(scriptDir, "bin/functions.R"))

library(tidyverse)

rawLenDis <- read_tsv(rawTcLenDisFile)

tcLenDis <- rawLenDis %>% filter(LD.type == "tcLenDis")

bgMinusLD <- subtractTcBG(tcLenDis, bgTime = normLDtime)

bgMinusReadSum <-
  bgMinusLD %>%
  select(-seqLen, -reads, -bg.reads, -bg.subtract) %>%
  distinct()

maxReads <-
  bgMinusReadSum %>%
  filter(time == maxTime, mir.type == "mature") %>%
  select(flybase_id, LD.type, read.sum) %>%
  rename(max.reads = read.sum)

bgMinusReadSum %>%
  unite(lendis, LD.type, timepoint, time, sep = ".") %>%
  spread(lendis, read.sum) %>%
  write_tsv(file.path(outDir, 'bgMinusTcReads.tsv'))

bgMinusReadSum %>%
  left_join(maxReads) %>%
  mutate(read.norm = read.sum / max.reads) %>% select(-read.sum, max.reads) %>%
  filter(time != normLDtime) %>%
  unite(lendis, LD.type, timepoint, time, sep = ".") %>%
  spread(lendis, read.norm) %>%
  write_tsv(file.path(outDir, 'bgMinusNormReads.tsv'))

bgMinusLD %>% write_tsv('bgMinusLenDis.tsv')

## Include separate mature and star PPM
matureLD <- convertLDtoWide(bgMinusLD, "mature")
starLD <- convertLDtoWide(bgMinusLD, "star")

matureLD %>% write_tsv('bgMinus.miR.lendis.tsv')
starLD %>% write_tsv('bgMinus.miRSTAR.lendis.tsv')
