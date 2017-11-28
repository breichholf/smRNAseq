setupRlibs <- function(R_lib){

  .libPaths( c( .libPaths(), R_lib ) )

  if (!require("Biostrings")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biostrings", suppressUpdates=TRUE)
  }

  if (!require("Rsamtools")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rsamtools", suppressUpdates=TRUE)
  }

  if (!require("BiocParallel")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("BiocParallel", suppressUpdates=TRUE)
  }

  if (!require("stringr")){
    install.packages("stringr", dependencies=TRUE, repos='http://cloud.r-project.org/')
  }

  if (!require("forcats")){
    install.packages("forcats", dependencies=TRUE, repos='http://cloud.r-project.org/')
  }

  if (!require("purrr")){
    install.packages("purrr", dependencies=TRUE, repos='http://cloud.r-project.org/')
  }

  if (!require("readr")){
    install.packages("readr", dependencies=TRUE, repos='http://cloud.r-project.org/')
  }

  if (!require("tibble")) {
    install.packages("tibble", dependencies=TRUE, repos='http://cloud.r-project.org/')
  }

  if (!require("tidyr")) {
    install.packages("tibble", dependencies=TRUE, repos='http://cloud.r-project.org/')
  }

  if (!require("cowplot")) {
    install.packages("cowplot", dependencies=TRUE, repos='http://cloud.r-project.org/')
  }

  if (!require("knitr")) {
    install.packages("knitr", dependencies=TRUE, repos='http://cloud.r-project.org/')
  }

  # if (!require("optparse")) {
  #   install.packages("optparse", dependencies=TRUE, repos='http://cloud.r-project.org/')
  #   library("optparse")
  # }

}

prepareFlybaseGFF <- function(gff, conversion, output) {
  require(tidyverse)
  require(stringr)

  # Read in file
  fbGFF <- read_tsv(opt$gff,
                    col_names = c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'dot', 'tags'),
                    col_types = list('c', 'c', 'c', 'n', 'n', 'c', 'c', 'c', 'c'))
  fbRows <- fbGFF %>% separate_rows(tags, sep = ";")
  fbIDs <-
    fbRows %>% filter(grepl("ID=FBtr", tags)) %>%
    mutate(id = str_sub(tags, 4, length(tags))) %>%
    select(type, id) %>% rename(parent.type = type)

  # miRNAs have "Parent=FBgn", so they're filtered out here
  fbGenicParents <-
    fbRows %>% filter(grepl("Parent=FBtr", tags)) %>%
    mutate(id = str_sub(tags, 8, length(tags))) %>%
    rename(source.type = type) %>% left_join(fbIDs, by = "id")

  fbGenic <-
    fbGenicParents %>% filter(!grepl("mito", chr) & grepl("exon|intron|UTR", source.type)) %>%
    mutate(print.start = ifelse(parent.type == "tRNA", ifelse(start < 20, 0, start - 20),
                                start),
           print.end = ifelse(parent.type == "tRNA",
                              end + 20,
                              end),
           print.name = ifelse(grepl("mRNA|pseudogene|ncRNA", parent.type),
                               paste(parent.type, source.type, sep = "_"),
                               parent.type),
           score = 0) %>%
    select(chr, print.start, print.end, print.name, score, strand)

  fbMirs <-
    fbGFF %>% filter(!grepl("mito", chr) & grepl("miRNA", type)) %>%
    mutate(score = 0) %>%
    rename(print.start = start, print.end = end, print.name = type) %>%
    select(chr, print.start, print.end, print.name, score, strand)

  # Concatenate genic and miRs
  concatBed <- rbind(fbGenic, fbMirs)

  # Convert Chr Names
  # Read in file first - first lines contain "#" which are comments.
  assemblyReport <- read_tsv(opt$nameconversion, comment = "#",
                             col_names = c('fb.name', 'sequence.role', 'assigned.molecule',
                                           'molecule.loaction.type', 'genbank.accn', 'relationship',
                                           'RefSeq.accn', 'Assembly.Unit', 'sequence.length', 'ucsc.name'),
                             col_types = list('c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'))

  # Select appropriate columns and remove leading "chr" from UCSC Names
  chrConversion <-
    assemblyReport %>%
    mutate(new.ucsc = str_replace(ucsc.name, "^chr", "")) %>%
    select(fb.name, new.ucsc) %>%
    mutate(fb.name = ifelse(new.ucsc == "M", "mitochondrion_genome", fb.name))

  # Convert Flybase Names to UCSC-compatible names
  outputBed <-
    concatBed %>% left_join(chrConversion, by = c('chr' = 'fb.name')) %>%
    select(new.ucsc, print.start, print.end, print.name, score, strand)

  # Write output file, will overwrite old files
  outputBed %>% write_tsv(opt$output, col_names = FALSE)
}

get.top.startpos <- function(id, align, sRNAreads, mirAnno = NULL, topn = 5, ...) {
  require(dplyr)
  require(stringr)
  require(Rsamtools)
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

  r.summary <-
    map.r %>%
    group_by(rname, pos, seqLen) %>% summarise(lenDis = n()) %>%
    group_by(rname, pos) %>% mutate(totalReads = sum(lenDis)) %>% ungroup() %>%
    mutate(flybase_id = as.character(rname))

  tc.summary <-
    map.r %>%
    dplyr::filter(!is.na(TC)) %>%
    group_by(rname, pos, seqLen) %>% summarise(tcLenDis = n()) %>%
    group_by(rname, pos) %>% mutate(tcReads = sum(tcLenDis)) %>% ungroup() %>%
    mutate(flybase_id = as.character(rname)) %>%
    dplyr::select(-rname)

  r.sum.pos <-
    r.summary %>%
    left_join(tc.summary, by = c("flybase_id", "pos", "seqLen")) %>%
    replace_na(list(totalReads = 0, lenDis = 0, tcReads = 0, tcLenDis = 0)) %>%
    left_join(mirAnno, by = "flybase_id") %>%
    dplyr::select(-rname, -loop)

  r.sum.arms <-
    r.sum.pos %>%
    dplyr::filter((pos >= `5p` - 5 & pos <= `5p` + 5) | (pos >= `3p` - 5 & pos <= `3p` + 5)) %>%
    mutate(arm.name = ifelse(pos >= `5p` - 5 & pos <= `5p` + 5,
                             paste0(str_sub(mir_name, 5, -1), "-5p"),
                             paste0(str_sub(mir_name, 5, -1), "-3p")))

  totalRName <- paste("totalReads", id, sep = ".")
  tcRName <- paste("tcReads", id, sep = ".")
  totalLDname <- paste("totalLenDis", id, sep = ".")
  tcLDname <- paste("tcLenDis", id, sep = ".")

  r.sum.max.pos <-
    r.sum.arms %>%
    group_by(arm.name) %>%
    top_n(n = topn, wt = totalReads) %>% ungroup() %>%
    mutate(totalReads = totalReads / sRNAreads * 1000000,
           tcReads = tcReads / sRNAreads * 1000000,
           lenDis = lenDis / sRNAreads * 1000000,
           tcLenDis = tcLenDis / sRNAreads * 1000000) %>%
    dplyr::rename_(.dots = setNames(c("totalReads", "tcReads", "lenDis", "tcLenDis"),
                                    c(totalRName, tcRName, totalLDname, tcLDname)))

  return(r.sum.max.pos)
}

convertToWide <- function(gatheredAllCounts, mirType) {
  require(dplyr)
  require(forcats)
  mirType <- enquo(mirType)

  output <-
    gatheredAllCounts %>%
    dplyr::filter(mir.type == !!mirType) %>%
    mutate(arm.name = fct_reorder(arm.name, desc(average.reads))) %>%
    dplyr::select(arm.name, pos, timepoint, reads) %>%
    spread(timepoint, reads)

  return(output)
}

mutsFromPileup <- function(flybase_id, pos, bamFile, timepoint, full.seq, minLen = 18, ...) {
  require(tibble)
  require(dplyr)
  require(stringr)
  require(Rsamtools)
  # Create tibble to merge in later on with reference sequence
  refSeq <- as_tibble(list("flybase_id" = flybase_id, "ref.seq" = full.seq))

  refSeqWpos <-
    refSeq %>%
    separate(ref.seq, paste("Pos", 1:str_length(full.seq), sep = "_"), sep = "\\B") %>%
    gather(pos, refNuc, matches("Pos_")) %>%
    separate(pos, c("pos", "idx"), sep = "_", convert = TRUE) %>%
    dplyr::select(-pos)

  # Get pileup
  pileupParams <- PileupParam(query_bins = seq(0,30), max_depth=10000000, min_mapq=0, min_base_quality=0)
  start.pos <- pos
  end.pos <- start.pos + 30
  id <- basename(bamFile) %>% str_sub(1, 5) # The first 5 characters of `bamFile` are the ID!
  scanParams <- ScanBamParam(flag = scanBamFlag(isMinusStrand = F),
                             which=GRanges(flybase_id, IRanges(start.pos, end.pos)))
  # Let's filter out N-containing reads for the entire file
  # We maybe could use this to filter for qualities <= 27, too?
  # That might "filter out any read with mutation quality <= 27", though!
  filterNs <- FilterRules(list(NoAmbigNucleotide = function(x) !grepl("N", x$seq)))
  filterBam <- filterBam(bamFile, tempfile(),
                         param = ScanBamParam(what = "seq",
                                              flag = scanBamFlag(isMinusStrand = F)),
                         filter = filterNs)

  # Actual command to get the pileup
  pileupResult <- pileup(filterBam, scanBamParam = scanParams, pileupParam = pileupParams)

  # Filter pileup to only get reads starting at our desired start position
  filteredRes <-
    pileupResult %>%
    dplyr::select(-which_label, -strand) %>%
    mutate(relPos = as.numeric(query_bin),
           flybase_id = as.character(seqnames), # Coerce factor to character to avoid warning later on
           timepoint = timepoint) %>%
    dplyr::select(-seqnames, -query_bin) %>%
    dplyr::filter(relPos == pos - min(pos) + 1, relPos <= minLen) %>%
    left_join(refSeqWpos, by = c("flybase_id", "pos" = "idx")) # Merge in `refSeqWpos` from above

  return(filteredRes)
}

onlyTCReads <- function(flybase_id, pos, bamFile, timepoint, full.seq, minLen = 18, ...) {
  require(tibble)
  require(dplyr)
  require(stringr)
  require(Rsamtools)
  # Create tibble to merge in later on with reference sequence
  refSeq <- as_tibble(list("flybase_id" = flybase_id, "ref.seq" = full.seq))

  refSeqWpos <-
    refSeq %>%
    separate(ref.seq, paste("Pos", 1:str_length(full.seq), sep = "_"), sep = "\\B") %>%
    gather(pos, refNuc, matches("Pos_")) %>%
    separate(pos, c("pos", "idx"), sep = "_", convert = TRUE) %>%
    dplyr::select(-pos)

  # Get pileup
  pileupParams <- PileupParam(query_bins = seq(0,30), max_depth=10000000, min_mapq=0, min_base_quality=0)
  start.pos <- pos
  end.pos <- start.pos + 30
  id <- basename(bamFile) %>% str_sub(1, 5) # The first 5 characters of `bamFile` are the ID!
  scanParams <- ScanBamParam(flag = scanBamFlag(isMinusStrand = F),
                             which=GRanges(flybase_id, IRanges(start.pos, end.pos)))
  # Let's filter out N-containing reads for the entire file
  # We maybe could use this to filter for qualities <= 27, too?
  # That might "filter out any read with mutation quality <= 27", though!
  onlyTCandNoNs <- FilterRules(list(NoAmbigNucleotide = function(x) !grepl("N", x$seq), TConly = function(x) x$TC == 1))
  filterBam <- filterBam(bamFile, tempfile(),
                         param = ScanBamParam(what = "seq",
                                              tag = c('TN', 'TC'),
                                              flag = scanBamFlag(isMinusStrand = F)),
                         filter = onlyTCandNoNs)
}

pileupParallelMuts <- function(groupedData, snow) {
  require(BiocParallel)
  require(dplyr)

  fbid <- groupedData$flybase_id
  bF <- unique(groupedData$bamFile)
  tp <- groupedData$timepoint
  pos <- groupedData$pos

  doOut <- bpmapply(doParallelPileup, miR = fbid, timepoint = tp, pos = pos,
                    MoreArgs = list(bamFile = bF, minLen = 18),
                    SIMPLIFY = FALSE, BPPARAM = snow)

  return(dplyr::bind_rows(doOut))
}

doParallelPileup <- function(miR, timepoint, pos, bamFile, minLen) {
  # This function will be called from dplyr do() in parallel using BiocParallel `bpmapply`
  # The function itself returns a cleaned data.frame of the pileup, which mapply wraps in a list
  # with one item for every bamFile.
  require(Rsamtools)
  require(dplyr)

  start.pos <- pos
  end.pos <- start.pos + 30

  pparam <- PileupParam(query_bins = seq(0,30), max_depth=10000000, min_mapq=0, min_base_quality=0)
  sparam <- ScanBamParam(flag = scanBamFlag(isMinusStrand = F),
                             which=GRanges(miR, IRanges(start.pos, end.pos)))

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
           timepoint = timepoint) %>%
    dplyr::select(-seqnames, -query_bin) %>%
    dplyr::filter(relPos == pos - min(pos) + 1, relPos <= minLen)

  return(filteredRes)
}

