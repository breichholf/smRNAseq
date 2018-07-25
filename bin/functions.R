setupRlibs <- function(R_lib){
  .libPaths( c( .libPaths(), R_lib ) )

  if(!require("pacman")) {
    install.packages("pacman", dependencies = TRUE, repos = 'http://cloud.r-project.org/')
  }

  p_load(rlang)

  p_install_version(
    c('rlang', 'tidyverse', 'cowplot', 'jsonlite', 'devtools', 'seqinr', 'sessioninfo',
      'biomaRt', 'Biostrings', 'Rsamtools', 'BiocParallel'),
    c('0.2.1', '1.2.1', '0.9.2', '1.5', '1.13.6', '3.4-5', '1.0.0',
      '2.32.1', '2.44.2', '1.28.0', '1.10.1')
  )

  p_load(sessioninfo)
}

getcfg <- function(json) {
  suppressMessages(require(jsonlite))
  suppressMessages(require(dplyr))
  suppressMessages(require(readr))

  cfg.info <- jsonlite::read_json(json)
  file.home <- cfg.info$base
  file.ncRNA <- cfg.info$ncBase
  mir.anno <- read_tsv(cfg.info$mir.anno)

  cfg.samples <- dplyr::bind_rows(cfg.info$samples)
  cfg.samples <- mutate(cfg.samples,
                        align = file.path(file.home, align),
                        ncBam = file.path(file.ncRNA, ncBam))

  return(list("samples" = cfg.samples, "anno" = mir.anno))
}

# 1) Reads bam
# 2) Filters out reads that don't map with in +/- 5 of annotated 5p and 3p arms
# 3) Counts all reads (and all TC reads with BQ>27) for given starting position
# 4) Assesses read lengths
# 5) Normalises reads to sRNAreads provided in cfg file
getAllCounts <- function(id, align, ncBam, sRNAreads, time, mirAnno = NULL, topn = 5, ...) {
  suppressMessages(require(tidyverse))
  suppressMessages(require(Rsamtools))
  mapInfo <- c("rname", "strand", "pos")
  mapParams <- ScanBamParam(what = c(mapInfo, "seq"), tag = c("TC", "TN"),
                            flag = scanBamFlag(isMinusStrand = FALSE, isUnmappedQuery = FALSE))
  filterNs <- FilterRules(list(NoAmbigNucleotide = function(x) !grepl("N", x$seq)))
  filterBam <- filterBam(align, tempfile(), filter = filterNs)
  bam <- scanBam(filterBam, param = mapParams)
  # Now this will ONLY handle files that have tags TC and TN, too!
  map.r <- dplyr::bind_cols(do.call(dplyr::bind_cols, bam[[1]][mapInfo]),
                            list("seqLen" = width(bam[[1]]$seq)),
                            do.call(dplyr::bind_cols, bam[[1]]$tag))

  # Sum up all reads with length X, to get total read count
  totalReadCounts <-
    map.r %>%
    group_by(rname, pos, seqLen) %>% summarise(lenDis = n()) %>%
    group_by(rname, pos) %>% mutate(totalReads = sum(lenDis)) %>% ungroup() %>%
    mutate(flybase_id = as.character(rname))

  # Sum up all reads where the custom TC flag was found
  tcReadCounts <-
    map.r %>%
    dplyr::filter(!is.na(TC)) %>%
    group_by(rname, pos, seqLen) %>% summarise(tcLenDis = n()) %>%
    group_by(rname, pos) %>% mutate(tcReads = sum(tcLenDis)) %>% ungroup() %>%
    mutate(flybase_id = as.character(rname)) %>%
    dplyr::select(-rname)

  # Join tc Read count in to total reads
  read.summary <-
    totalReadCounts %>%
    left_join(tcReadCounts, by = c("flybase_id", "pos", "seqLen")) %>%
    replace_na(list(totalReads = 0, lenDis = 0, tcReads = 0, tcLenDis = 0)) %>%
    left_join(mirAnno, by = "flybase_id") %>%
    dplyr::select(-rname)

  # Only keep reads that are within +/-10nt of the suggested 5p/3p arm positions
  read.summary.closestArms <-
    read.summary %>%
    dplyr::filter((pos >= `5p` - 10 & pos <= `5p` + 10) | (pos >= `3p` - 10 & pos <= `3p` + 10)) %>%
    mutate(arm.name = ifelse(pos >= `5p` - 10 & pos <= `5p` + 10,
                             paste0(str_sub(mir_name, 5, -1), "-5p"),
                             paste0(str_sub(mir_name, 5, -1), "-3p")))

  totalRName <- paste("totalReads", id, time, sep = ".")
  tcRName <- paste("tcReads", id, time, sep = ".")
  totalLDname <- paste("totalLenDis", id, time, sep = ".")
  tcLDname <- paste("tcLenDis", id, time, sep = ".")

  # Keep the `topn` (default: 5) most frequently used start positions
  read.summary.topnPos <-
    read.summary.closestArms %>%
    group_by(arm.name) %>%
    top_n(n = topn, wt = totalReads) %>% ungroup() %>%
    mutate(totalReads = totalReads / sRNAreads * 1000000,
           tcReads = tcReads / sRNAreads * 1000000,
           lenDis = lenDis / sRNAreads * 1000000,
           tcLenDis = tcLenDis / sRNAreads * 1000000) %>%
    dplyr::rename_(.dots = setNames(c("totalReads", "tcReads", "lenDis", "tcLenDis"),
                                    c(totalRName, tcRName, totalLDname, tcLDname)))

  return(read.summary.topnPos)
}

pileupParallelMuts <- function(groupedData, mc.param, minLen) {
  suppressMessages(require(BiocParallel))
  suppressMessages(require(dplyr))

  fbid <- groupedData$flybase_id
  bF <- unique(groupedData$bamFile)
  tp <- groupedData$timepoint
  t <- groupedData$time
  pos <- groupedData$pos
  mir.type <- groupedData$mir.type

  doOut <- bpmapply(doParallelPileup, miR = fbid, timepoint = tp, time = t, pos = pos,
                    mir.type = mir.type, MoreArgs = list(bamFile = bF, minLen = minLen),
                    SIMPLIFY = FALSE, BPPARAM = mc.param)

  return(dplyr::bind_rows(doOut))
}

doParallelPileup <- function(miR, timepoint, time, pos, mir.type, bamFile, minLen) {
  # This function will be called from dplyr do() in parallel using BiocParallel `bpmapply`
  # The function itself returns a cleaned data.frame of the pileup, which mapply wraps in a list
  # with one item for every bamFile.
  suppressMessages(require(Rsamtools))
  suppressMessages(require(dplyr))

  start.pos <- pos
  end.pos <- start.pos + 30

  pparam <- PileupParam(query_bins = seq(0,30), max_depth = 50000000, min_mapq = 0, min_base_quality = 0)
  sparam <- ScanBamParam(flag = scanBamFlag(isMinusStrand = F),
                         which = GRanges(miR, IRanges(start.pos, end.pos)))

  filterNs <- FilterRules(list(NoAmbigNucleotide = function(x) !grepl("N", x$seq)))
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
           time = time,
           mir.type = mir.type,
           start.pos = start.pos) %>%
    dplyr::select(-seqnames, -query_bin) %>%
    dplyr::filter(relPos == pos - min(pos) + 1, relPos <= minLen)

  return(filteredRes)
}

# 1) Reads bam
# 3) Counts all reads (and all TC reads with BQ>27) for given starting position
# 4) Assesses read lengths
# 5) Normalises reads to sRNAreads provided in cfg file
getNcRNACounts <- function(id, align, ncBam, sRNAreads, time, topn = 10, minLen = 18, maxLen = 30, ...) {
  suppressMessages(require(tidyverse))
  suppressMessages(require(Rsamtools))
  mapInfo <- c("rname", "strand", "pos")
  mapParams <- ScanBamParam(what = c(mapInfo, "seq"), tag = c("TC", "TN"),
                            flag = scanBamFlag(isMinusStrand = FALSE, isUnmappedQuery = FALSE))
  filterNs <- FilterRules(list(NoAmbigNucleotide = function(x) !grepl("N", x$seq),
                               MinLen = function(x) nchar(x$seq) >= minLen,
                               MaxLen = function(x) nchar(x$seq) <= maxLen))
  filterBam <- filterBam(ncBam, tempfile(), filter = filterNs)
  bam <- scanBam(filterBam, param = mapParams)
  # Now this will ONLY handle files that have tags TC and TN, too!
  map.r <- dplyr::bind_cols(do.call(dplyr::bind_cols, bam[[1]][mapInfo]),
                            list("seqLen" = width(bam[[1]]$seq)),
                            do.call(dplyr::bind_cols, bam[[1]]$tag))

  # Sum up all reads with length X, to get total read count
  totalReadCounts <-
    map.r %>%
    group_by(rname, pos, seqLen) %>% summarise(lenDis = n()) %>%
    group_by(rname, pos) %>% mutate(totalReads = sum(lenDis)) %>% ungroup() %>%
    mutate(locus = as.character(rname))

  # Sum up all reads where the custom TC flag was found
  tcReadCounts <-
    map.r %>%
    dplyr::filter(!is.na(TC)) %>%
    group_by(rname, pos, seqLen) %>% summarise(tcLenDis = n()) %>%
    group_by(rname, pos) %>% mutate(tcReads = sum(tcLenDis)) %>% ungroup() %>%
    mutate(locus = as.character(rname)) %>%
    dplyr::select(-rname)

  # Join tc Read count in to total reads
  read.summary <-
    totalReadCounts %>%
    left_join(tcReadCounts, by = c("locus", "pos", "seqLen")) %>%
    replace_na(list(totalReads = 0, lenDis = 0, tcReads = 0, tcLenDis = 0)) %>%
    dplyr::select(-rname)

  totalRName <- paste("totalReads", id, time, sep = ".")
  tcRName <- paste("tcReads", id, time, sep = ".")
  totalLDname <- paste("totalLenDis", id, time, sep = ".")
  tcLDname <- paste("tcLenDis", id, time, sep = ".")

  # Keep the `topn` (default: 10) most frequently used start positions for each locus
  read.summary.topnPos <-
    read.summary %>%
    group_by(locus) %>%
    top_n(n = topn, wt = totalReads) %>% ungroup() %>%
    mutate(totalReads = totalReads / sRNAreads * 1000000,
           tcReads = tcReads / sRNAreads * 1000000,
           lenDis = lenDis / sRNAreads * 1000000,
           tcLenDis = tcLenDis / sRNAreads * 1000000) %>%
    dplyr::rename_(.dots = setNames(c("totalReads", "tcReads", "lenDis", "tcLenDis"),
                                    c(totalRName, tcRName, totalLDname, tcLDname)))

  return(read.summary.topnPos)
}

pileupMutsLenRestrict <- function(groupedData, minLen, maxLen, mc.param) {
  suppressMessages(require(BiocParallel))
  suppressMessages(require(dplyr))

  locus <- groupedData$locus
  bF <- unique(groupedData$bamFile)
  tp <- groupedData$timepoint
  t <- groupedData$time
  pos <- groupedData$pos

  doOut <- bpmapply(parallelMutsLenRestrict, locus = locus, timepoint = tp, time = t, pos = pos,
                    MoreArgs = list(bamFile = bF, minLen = minLen, maxLen = maxLen),
                    SIMPLIFY = FALSE, BPPARAM = mc.param)

  return(dplyr::bind_rows(doOut))
}

parallelMutsLenRestrict <- function(locus, timepoint, time, pos, biotype, bamFile, minLen, maxLen) {
  suppressMessages(require(Rsamtools))
  suppressMessages(require(dplyr))

  start.pos <- pos
  end.pos <- start.pos + 30

  if (maxLen < 30) { lastBin <- maxLen }
  else { lastBin <- 30 }

  pparam <- PileupParam(query_bins = seq(0,lastBin), max_depth = 50000000, min_mapq = 0, min_base_quality = 0)
  sparam <- ScanBamParam(flag = scanBamFlag(isMinusStrand = F),
                         which = GRanges(locus, IRanges(start.pos, end.pos)))

  filterNsAndLens <- FilterRules(list(NoAmbigNucleotide = function(x) !grepl("N", x$seq),
                                      MinLen = function(x) nchar(x$seq) >= minLen,
                                      MaxLen = function(x) nchar(x$seq) <= maxLen))

  filterBam <- filterBam(bamFile, tempfile(),
                         param = ScanBamParam(what = "seq",
                                              flag = scanBamFlag(isMinusStrand = F)),
                         filter = filterNsAndLens)

  pileupResult <- pileup(filterBam, scanBamParam = sparam, pileupParam = pparam)

  filteredRes <-
    pileupResult %>%
    dplyr::select(-which_label, -strand) %>%
    mutate(relPos = as.numeric(query_bin),
           locus = as.character(seqnames), # Coerce factor to character to avoid warning later on
           timepoint = timepoint,
           time = time,
           start.pos = start.pos) %>%
    dplyr::select(-seqnames, -query_bin) %>%
    dplyr::filter(relPos == pos - min(pos) + 1, relPos <= maxLen) # Retain all mutations

  return(filteredRes)
}
