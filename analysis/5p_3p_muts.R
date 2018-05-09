#!/usr/bin/env Rscript

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
scriptDir <- as.character(args[1])
R_libs <- as.character(args[2])
nCores <- as.numeric(args[3])
jsonFile <- as.character(args[4])
ldFile <- as.character(args[5])
preMirFastaFile <- as.character(args[6])

source(file.path(scriptDir, 'bin/functions.R'))

setupRlibs(R_libs)

pileupAllMuts <- function(groupedData, mc.param) {
  suppressMessages(require(BiocParallel))
  suppressMessages(require(dplyr))
  
  fbid <- groupedData$flybase_id
  bF <- unique(groupedData$bamFile)
  tp <- groupedData$timepoint
  pos <- groupedData$pos
  mL <- groupedData$seqLen
  
  doOut <- bpmapply(doAllMuts, miR = fbid, timepoint = tp, pos = pos, maxLen = mL,
                    MoreArgs = list(bamFile = bF), SIMPLIFY = FALSE,
                    BPPARAM = mc.param)
  
  return(dplyr::bind_rows(doOut))
}

doAllMuts <- function(miR, timepoint, pos, maxLen, bamFile) {
  # This function will be called from dplyr do() in parallel using BiocParallel `bpmapply`
  # The function itself returns a cleaned data.frame of the pileup, which mapply wraps in a list
  # with one item for every bamFile.
  suppressMessages(require(Rsamtools))
  suppressMessages(require(dplyr))
  
  start.pos <- pos
  end.pos <- start.pos + 30
  
  pparam <- PileupParam(query_bins = seq(0,30), max_depth=50000000, min_mapq=0, min_base_quality=0)
  sparam <- ScanBamParam(flag = scanBamFlag(isMinusStrand = F),
                         which=GRanges(miR, IRanges(start.pos, end.pos)))
  
  filterNsMaxLen <- FilterRules(list(NoAmbigNucleotide = function(x) !grepl("N", x$seq),
                                     SingleLen = function(x) nchar(x$seq) == maxLen))
  filterBam <- filterBam(bamFile, tempfile(),
                         param = ScanBamParam(what = "seq",
                                              flag = scanBamFlag(isMinusStrand = F)),
                         filter = filterNs)
  
  pileupResult <- pileup(filterBam, scanBamParam = sparam, pileupParam = pparam)
  
  filteredRes <-
    pileupResult %>%
    dplyr::select(-which_label, -strand) %>%
    mutate(relPos = as.numeric(query_bin),
           flybase_id = as.character(seqnames), # Coerce factor to character to avoid warning later on
           timepoint = timepoint,
           start.pos = start.pos) %>%
    dplyr::select(-seqnames, -query_bin) %>%
    dplyr::filter(relPos == pos - min(pos) + 1, relPos <= maxLen)
  
  return(filteredRes)
}

library(tidyverse)
library(Biostrings)
library(BiocParallel)

sessionInfo()

cfg <- getcfg(jsonFile)

mirPositions <- read_tsv(ldFile)

preMirFasta <- readDNAStringSet(preMirFastaFile)
preMirTbl <- as_tibble(list('flybase_id' = names(preMirFasta), 'full.seq' = paste(preMirFasta)))

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
  left_join(preMirTbl)

topPosWlens <-
  mirPosWFiles %>%
  filter(LD.type == "totalLenDis") %>%
  gropu_by(arm.name, pos, timepoint) %>%
  top_n(1, reads) %>%
  gropu_by(flybase_id, arm.name, pos, seqLen) %>%
  summarise(counter = n()) %>%
  group_by(arm.name, pos) %>%
  filter(counter == max(counter), counter > 1)

distinctFiles <-
  mirPosWFiles %>%
  select(flybase_id, arm.name, pos, bamFile, timepoint) %>%
  distinct()

doOut <- bpmapply(doAllMuts, miR = fbid, timepoint = tp, pos = pos, maxLen = mL,
                  MoreArgs = list(bamFile = bF), SIMPLIFY = FALSE,
                  BPPARAM = mc.param)

# This is a bit hackey, as on our cluster, the nCores provided is half what R can use by hyperthreading
nProcs <- nCores * 2
mc.param <- MulticoreParam(workers = nProcs, type = 'FORK')

mirsWmuts <-
  topPosWlens %>%
  left_join(distinctFiles) %>%
  group_by(bamFile) %>%
  do(muts = pileupAllMuts(groupedData = ., mc.param = mc.param)) %>%
  unnest(muts) %>%
  dplyr::select(-bamFile) %>%
  left_join(tidyRefNucs, by = c('flybase_id', 'pos' = 'idx')) %>%
  left_join(mirPosWFiles %>% dplyr::select(-bamFile), by = c('flybase_id', 'timepoint', 'time', 'mir.type', 'start.pos' = 'pos'))

mirsWmuts %>% write_tsv('miRs.entireLenMuts.tsv')