library(Biostrings)
library(tidyverse)
library(jsonlite)

pMFF <- '/Users/reichholf/bioinfo/data/vh/20171020_S2_LNA-transfection/results/bowtie-sorted/ext_hairpins/hairpin.fa'
preMirFasta <- readDNAStringSet(pMFF)
preMirTbl <- as_tibble(list("flybase_id" = names(preMirFasta), "full.seq" = paste(preMirFasta)))

mutsFromPileup <- function(flybase_id, arm.name, pos, bamFile, timepoint, full.seq, minLen = 18, ...) {
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

mutsToCode <- function(mutDF, mirMetaInfoDF) {
  mutCodes <-
    mutDF %>%
    spread(nucleotide, count) %>%
    replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
    gather(nucleotide, count, A:T) %>%
    mutate(mutCode = ifelse(refNuc != nucleotide,
                            paste(refNuc, nucleotide, sep = ">"),
                            refNuc)) %>%
    group_by(flybase_id, timepoint, pos) %>%
    mutate(depth = sum(count), mutFract = count / depth) %>%
    dplyr::filter(grepl(">", mutCode)) %>%
    arrange(flybase_id, timepoint, relPos, mutCode) %>%
    ungroup() %>%
    mutate(first.pos = pos - relPos + 1) %>%
    left_join(mirMetaInfoDF)
  
  return(mutCodes)
}

makeMutsWide <- function(mutCodeDF) {
  wideMuts <-
    mutCodeDF %>%
    group_by(arm.name, timepoint, mutCode) %>%
    mutate(relMutCount = rank(relPos)) %>%
    mutate(relMut = paste(mutCode, relMutCount, sep = "_"),
           relMut = str_replace(relMut, ">", "")) %>% ungroup() %>%
    dplyr::select(arm.name, first.pos, seed, UCount, reads, average.reads, mir.type, timepoint, mutFract, relMut) %>%
    spread(relMut, mutFract) %>%
    arrange(mir.type, desc(average.reads), timepoint)
  
  return(wideMuts)
}

inFold <- '/Volumes/ameres/Reichholf/sequencing/vh/20171020_S2_LNA-transfection/results/stats/'
jsonFile <- '/Volumes/ameres/Reichholf/sequencing/vh/20171020_S2_LNA-transfection/results/json/samples.json'
readFile <- '/Volumes/ameres/Reichholf/sequencing/vh/20171020_S2_LNA-transfection/20171024_LNA-transfection_smallRNAreads.txt'
readOverview <- read_tsv(readFile)

selectMirs <- c("dme-bantam", "dme-mir-184")
scrLibs <- c(57518, 57519, 57520, 57521)
banLibs <- c(57522, 57523, 57524, 57525)
mir184Libs <- c(57526, 57527, 57528, 57529)

ctrlMirs <- c("dme-mir-14", "dme-mir-34", "dme-mir-998", "dme-mir-988")

bed <- read_tsv('/Users/reichholf/bioinfo/data/vh/20171020_S2_LNA-transfection/results/bowtie-sorted/ext_hairpins/hairpin_plus20nt.bed',
                col_names = c("chr", "start", "end", 'flybase_id', "score", "strand"))

cfg.info <- jsonlite::read_json(jsonFile)
file.home <- cfg.info$base
mir.anno <- read_tsv('/Volumes/ameres/Reichholf/sequencing/ref-seq/r6.08/dmel_r6.08_FB2015_05/fasta/mir_arm_positions-r6.08.tsv')

cfg.samples <- dplyr::bind_rows(cfg.info$samples)
cfg.samples <- mutate(cfg.samples, align = file.path(file.home, align))

tPCF <- file.path(inFold, 'topPositionCounts.tsv')
aCF <- file.path(inFold, 'allCounts.tsv')

tPC <- read_tsv(tPCF)
aC <- read_tsv(aCF)

cfg.local <-
  cfg.samples %>%
  mutate(align = str_replace(align, "/scratch/brian/lna/", '/Users/reichholf/bioinfo/data/vh/20171020_S2_LNA-transfection/'))

mirBodyLength <- 18

mirSubset <-
  tPC %>%
  filter(mir_name %in% selectMirs) %>%
  arrange(mir_name, desc(average.reads), timepoint)

mirSubsetFull <-
  mirSubset %>%
  left_join(cfg.local, by = c('timepoint' = 'id')) %>%
  mutate(bamFile = file.path(align)) %>%
  left_join(preMirTbl)

mirSubsetRedux <-
  mirSubsetFull %>%
  select(pos, flybase_id, mir_name, arm.name, timepoint, reads, average.reads, mir.type, seed, UCount, sRNAreads) %>%
  rename(first.pos = pos)

mirScrSubset <- mirSubsetFull %>% filter(timepoint %in% scrLibs)
mirBanSubset <- mirSubsetFull %>% filter(timepoint %in% banLibs)
mir184Subset <- mirSubsetFull %>% filter(timepoint %in% mir184Libs)

mirScrMuts <- mirScrSubset %>% pmap_dfr(mutsFromPileup, minLen = mirBodyLength)
mirBanMuts <- mirBanSubset %>% pmap_dfr(mutsFromPileup, minLen = mirBodyLength)
mir184Muts <- mir184Subset %>% pmap_dfr(mutsFromPileup, minLen = mirBodyLength)

mirScrCodes <- mutsToCode(mirScrMuts, mirSubsetRedux)
mirBanCodes <- mutsToCode(mirBanMuts, mirSubsetRedux)
mir184Codes <- mutsToCode(mir184Muts, mirSubsetRedux)

mirScrWide <- makeMutsWide(mirScrCodes)
mirBanWide <- makeMutsWide(mirBanCodes)
mir183Wide <- makeMutsWide(mir184Codes)

mirScrCodes <-
  mirScrMuts %>%
  spread(nucleotide, count) %>%
  replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
  gather(nucleotide, count, A:T) %>%
  mutate(mutCode = ifelse(refNuc != nucleotide,
                          paste(refNuc, nucleotide, sep = ">"),
                          refNuc)) %>%
  group_by(flybase_id, timepoint, pos) %>%
  mutate(depth = sum(count), mutFract = count / depth) %>%
  dplyr::filter(grepl(">", mutCode)) %>%
  arrange(flybase_id, timepoint, relPos, mutCode) %>%
  ungroup() %>%
  mutate(first.pos = pos - relPos + 1) %>%
  left_join(mirSubsetRedux)

mirScrWide <-
  mirScrCodes %>%
  group_by(arm.name, timepoint, mutCode) %>%
  mutate(relMutCount = rank(relPos)) %>%
  mutate(relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", "")) %>% ungroup() %>%
  dplyr::select(arm.name, first.pos, seed, UCount, reads, average.reads, mir.type, timepoint, mutFract, relMut) %>%
  spread(relMut, mutFract) %>%
  arrange(mir.type, desc(average.reads), timepoint)

mirScrWide %>% write_tsv(file.path(inFold, 'scr-LNA_ban-mir-184_muts.tsv'))
mirBanWide %>% write_tsv(file.path(inFold, 'ban-LNA_ban-mir-184_muts.tsv'))
mir184Wide %>% write_tsv(file.path(inFold, 'mir184-LNA_ban-mir-184_muts.tsv'))

mirCtrlSubset <-
  tPC %>%
  filter(mir_name %in% ctrlMirs) %>%
  arrange(mir_name, desc(average.reads), timepoint) %>%
  left_join(cfg.local, by = c('timepoint' = 'id')) %>%
  mutate(bamFile = file.path(align)) %>%
  left_join(preMirTbl)

mirCtrlRedux <-
  mirCtrlSubset %>%
  select(pos, flybase_id, mir_name, arm.name, timepoint, reads, average.reads, mir.type, seed, UCount, sRNAreads) %>%
  rename(first.pos = pos)

# mirCtrlScrCodes <-
#   mirCtrlScrMuts %>%
#   spread(nucleotide, count) %>%
#   replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
#   gather(nucleotide, count, A:T) %>%
#   mutate(mutCode = ifelse(refNuc != nucleotide,
#                           paste(refNuc, nucleotide, sep = ">"),
#                           refNuc)) %>%
#   group_by(flybase_id, timepoint, pos) %>%
#   mutate(depth = sum(count), mutFract = count / depth) %>%
#   dplyr::filter(grepl(">", mutCode)) %>%
#   arrange(flybase_id, timepoint, relPos, mutCode) %>%
#   ungroup() %>%
#   mutate(first.pos = pos - relPos + 1) %>%
#   left_join(mirCtrlRedux)
# 
# mirCtrlScrWide <-
#   mirCtrlScrCodes %>%
#   group_by(arm.name, timepoint, mutCode) %>%
#   mutate(relMutCount = rank(relPos)) %>%
#   mutate(relMut = paste(mutCode, relMutCount, sep = "_"),
#          relMut = str_replace(relMut, ">", "")) %>% ungroup() %>%
#   dplyr::select(arm.name, first.pos, seed, UCount, reads, average.reads, mir.type, timepoint, mutFract, relMut) %>%
#   spread(relMut, mutFract) %>%
#   arrange(mir.type, desc(average.reads), timepoint)
# 
# mirCtrlBanCodes <-
#   mirCtrlBanMuts %>%
#   spread(nucleotide, count) %>%
#   replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
#   gather(nucleotide, count, A:T) %>%
#   mutate(mutCode = ifelse(refNuc != nucleotide,
#                           paste(refNuc, nucleotide, sep = ">"),
#                           refNuc)) %>%
#   group_by(flybase_id, timepoint, pos) %>%
#   mutate(depth = sum(count), mutFract = count / depth) %>%
#   dplyr::filter(grepl(">", mutCode)) %>%
#   arrange(flybase_id, timepoint, relPos, mutCode) %>%
#   ungroup() %>%
#   mutate(first.pos = pos - relPos + 1) %>%
#   left_join(mirCtrlRedux)
# 
# mirCtrlBanWide <-
#   mirCtrlBanCodes %>%
#   group_by(arm.name, timepoint, mutCode) %>%
#   mutate(relMutCount = rank(relPos)) %>%
#   mutate(relMut = paste(mutCode, relMutCount, sep = "_"),
#          relMut = str_replace(relMut, ">", "")) %>% ungroup() %>%
#   dplyr::select(arm.name, first.pos, seed, UCount, reads, average.reads, mir.type, timepoint, mutFract, relMut) %>%
#   spread(relMut, mutFract) %>%
#   arrange(mir.type, desc(average.reads), timepoint)
# 
# mirCtrl184Codes <-
#   mirCtrl184Muts %>%
#   spread(nucleotide, count) %>%
#   replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
#   gather(nucleotide, count, A:T) %>%
#   mutate(mutCode = ifelse(refNuc != nucleotide,
#                           paste(refNuc, nucleotide, sep = ">"),
#                           refNuc)) %>%
#   group_by(flybase_id, timepoint, pos) %>%
#   mutate(depth = sum(count), mutFract = count / depth) %>%
#   dplyr::filter(grepl(">", mutCode)) %>%
#   arrange(flybase_id, timepoint, relPos, mutCode) %>%
#   ungroup() %>%
#   mutate(first.pos = pos - relPos + 1) %>%
#   left_join(mirCtrlRedux)
# 
# mirCtrl184Wide <-
#   mirCtrl184Codes %>%
#   group_by(arm.name, timepoint, mutCode) %>%
#   mutate(relMutCount = rank(relPos)) %>%
#   mutate(relMut = paste(mutCode, relMutCount, sep = "_"),
#          relMut = str_replace(relMut, ">", "")) %>% ungroup() %>%
#   dplyr::select(arm.name, first.pos, seed, UCount, reads, average.reads, mir.type, timepoint, mutFract, relMut) %>%
#   spread(relMut, mutFract) %>%
#   arrange(mir.type, desc(average.reads), timepoint)

mirCtrlScr <- filter(mirCtrlSubset, timepoint %in% scrLibs)
mirCtrlBan <- filter(mirCtrlSubset, timepoint %in% banLibs)
mirCtrl184 <- filter(mirCtrlSubset, timepoint %in% mir184Libs)

mirCtrlScrMuts <- mirCtrlScr %>% pmap_dfr(mutsFromPileup, minLen = mirBodyLength)
mirCtrlBanMuts <- mirCtrlBan %>% pmap_dfr(mutsFromPileup, minLen = mirBodyLength)
mirCtrl184Muts <- mirCtrl184 %>% pmap_dfr(mutsFromPileup, minLen = mirBodyLength)

mirCtrlScrCodes <- mutsToCode(mirCtrlScrMuts, mirCtrlRedux)
mirCtrlBanCodes <- mutsToCode(mirCtrlBanMuts, mirCtrlRedux)
mirCtrl184Codes <- mutsToCode(mirCtrl184Muts, mirCtrlRedux)

mirCtrlScrWide <- makeMutsWide(mirCtrlScrCodes)
mirCtrlBanWide <- makeMutsWide(mirCtrlBanCodes)
mirCtrl184Wide <- makeMutsWide(mirCtrl184Codes)

mirCtrlScrWide %>% write_tsv(file.path(inFold, "scr-LNA_ctrl-miR_muts.tsv"))
mirCtrlBanWide %>% write_tsv(file.path(inFold, "ban-LNA_ctrl-miR_muts.tsv"))
mirCtrl184Wide %>% write_tsv(file.path(inFold, "mir184-LNA_ctrl-miR_muts.tsv"))

mirAllRedux <- dplyr::bind_rows(mirSubsetRedux, mirCtrlRedux)

mirAllScrMuts <- dplyr::bind_rows(mirScrMuts, mirCtrlScrMuts)
mirAllBanMuts <- dplyr::bind_rows(mirBanMuts, mirCtrlBanMuts)
mirAll184Muts <- dplyr::bind_rows(mir184Muts, mirCtrl184Muts)

mirAllScrCodes <- mutsToCode(mirAllScrMuts, mirAllRedux)
mirAllBanCodes <- mutsToCode(mirAllBanMuts, mirAllRedux)
mirAll184Codes <- mutsToCode(mirAll184Muts, mirAllRedux)

mirAllScrWide <- makeMutsWide(mirAllScrCodes)
mirAllBanWide <- makeMutsWide(mirAllBanCodes)
mirAll184Wide <- makeMutsWide(mirAll184Codes)

mirAllScrWide %>% write_tsv(file.path(inFold, "scr-LNA_all-miR_muts.tsv"))
mirAllBanWide %>% write_tsv(file.path(inFold, "ban-LNA_all-miR_muts.tsv"))
mirAll184Wide %>% write_tsv(file.path(inFold, "mir184-LNA_all-miR_muts.tsv"))