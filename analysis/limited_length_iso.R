ld_file <- '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/R6038 - QIAseq spike-in/raw/filteredLenDis.tsv'

ld_df <- read_tsv(ld_file)

bam <- '~/bioinfo/bamtest/65081_TACAGCGTTCAGAG_wt_nos4U.trimmed.TCtagged_hairpin.sorted.bam'

start.pos <- 52
end.pos <- start.pos + 30

mirLen <- 23

miR <- 'FBgn0262451'

pparam <- PileupParam(query_bins = seq(0,mirLen), max_depth=50000000, min_mapq=0, min_base_quality=0)

sparam <- ScanBamParam(flag = scanBamFlag(isMinusStrand = F),
                       which=GRanges(miR, IRanges(start.pos, end.pos)))

filterNs <- FilterRules(list(NoAmbigNucleotide = function(x) !grepl("N", x$seq),
                             SingleLen = function(x) nchar(x$seq) == mirLen))


filterBam <- filterBam(bam, tempfile(),
                       param = ScanBamParam(what = "seq",
                                            flag = scanBamFlag(isMinusStrand = F)),
                       filter = filterNs)

pileupResult_cbin <- pileup(filterBam, scanBamParam = sparam, pileupParam = pparam)

pileupResult_filterLen %>%
  dplyr::select(-which_label, -strand) %>%
  mutate(relPos = as.numeric(query_bin),
         flybase_id = as.character(seqnames),
         start.pos = start.pos) %>%
  dplyr::select(-seqnames, -query_bin) %>%
  dplyr::filter(relPos == pos - min(pos) + 1) %>%
  group_by(relPos) %>%
  mutate(depth = sum(count)) %>%
  filter(relPos >= 20)

  