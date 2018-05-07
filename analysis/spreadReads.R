#!/usr/bin/env Rscript

# Command line arguments
args = commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
rawLDFile <- as.character(args[2])
topPosFile <- as.character(args[3])
bgTime <- as.numeric(args[4])
normTime <- as.numeric(args[5])
outDir <- as.character(args[6])

source(file.path(scriptDir, "bin/functions.R"))

library(tidyverse)

rawLenDis <- read_tsv(rawLDFile)
topPositions <- read_tsv(topPosFile)

topPosLenDis <-
  rawLenDis %>%
  left_join(topPositions %>%
    select(pos, flybase_id, average.reads) %>%
    distinct()) %>%
  filter(average.reads >= 50)

steadyStateLD <-
  topPosLenDis %>%
  filter(LD.type == "totalLenDis") %>%
  mutate(LD.type = "stdyStateReads")

tcLenDis <-
  topPosLenDis %>%
  filter(LD.type == "tcLenDis")

bgMinusLD <-
  subtractTcBG(tcLenDis, bgTime = bgTime) %>%
  mutate(LD.type = "tcReads")

bgMinusReadSum <-
  bgMinusLD %>%
  select(-seqLen, -reads, -mir_name, -`5p`, -`3p`) %>%
  distinct()

normReads <-
  bgMinusLD %>%
  filter(time == normTime, mir.type == "mature") %>%
  select(flybase_id, LD.type, read.sum) %>% distinct() %>%
  dplyr::rename(norm.reads = read.sum)

####
# Output section
####
#### Steady state length distributions -- Raw
steadyStateLD %>% write_tsv(file.path(outDir, 'steadyState.raw.lendis.tsv'))

#    converted to wide and split in to separate output for mature and star
#    FOR NIBBLER PLOTS
convertLDtoWide(steadyStateLD, "mature") %>% write_tsv(file.path(outDir, 'steadyState.miR.lendis.tsv'))
convertLDtoWide(steadyStateLD, "star") %>% write_tsv(file.path(outDir, 'steadyState.miRSTAR.lendis.tsv'))

#### TC Reads
#    Background subtracted total reads
#    Tidy format
#    FOR PLOTS, calculation of kBIO FITS and REPLICATES
bgMinusReadSum %>% write_tsv(file.path(outDir, 'bgMinus.tidy.raw.TcReads.tsv'))

#    FOR SCISSOR PLOTS global
bgMinusScissor <-
  bgMinusReadSum %>%
  unite(lendis, LD.type, timepoint, time, sep = ".") %>%
  spread(lendis, read.sum) %>%
  arrange(mir.type, desc(average.reads))

bgMinusScissor %>% write_tsv(file.path(outDir, 'bgMinus.raw.TcReads.tsv'))

bgMinusScissor %>%
  filter(mir.type == "mature") %>%
  write_tsv(file.path(outDir, 'bgMinus.miR.TcReads.tsv'))

bgMinusScissor %>%
  filter(mir.type == "star") %>%
  write_tsv(file.path(outDir, 'bgMinus.miRSTAR.TcReads.tsv'))

#    set mature miR Read count at normTime to 1
#    convert to spread
#    FOR SCISSOR PLOTS comparison amongst replicates
bgMinusReadSum %>%
  left_join(normReads) %>% select(-flybase_id) %>%
  mutate(read.norm = read.sum / norm.reads) %>% select(-read.sum, -norm.reads) %>%
  filter(time != bgTime) %>%
  unite(lendis, LD.type, timepoint, time, sep = ".") %>%
  spread(lendis, read.norm) %>%
  arrange(mir.type, desc(average.reads)) %>%
  write_tsv(file.path(outDir, 'bgMinus.TcReads.Normalised.tsv'))

#    Raw Background subtracted length distribution data
bgMinusLD %>%
  select(-flybase_id, -mir_name, -`5p`, -`3p`) %>%
  arrange(mir.type, desc(average.reads), timepoint, seqLen) %>%
  write_tsv(file.path(outDir, 'bgMinus.raw.TC-lendis.tsv'))

#   Separate mature and star length distributions
#   miR and miRSTAR length distributions split
#   FOR NIBBLER PLOTS
convertTcLDtoWide(bgMinusLD, "mature") %>% write_tsv(file.path(outDir, 'bgMinus.miR.TC-lendis.tsv'))
convertTcLDtoWide(bgMinusLD, "star") %>% write_tsv(file.path(outDir, 'bgMinus.miRSTAR.TC-lendis.tsv'))
