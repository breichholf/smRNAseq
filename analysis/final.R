library(tidyverse)
library(broom)
library(minpack.lm)
library(nlsr)

baseFolder <- "/Volumes/ameres/Reichholf/sequencing/"
baseOutput <- "~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/final_analysis_2018-11/"

rep1_folder <- paste0(baseFolder, "20161102_Ago2KO_slam-pulse/single_1811_nfrerun/stats")
rep2_folder <- paste0(baseFolder, "20161205_Ago1KO_slam-pulse_R2/single_1811_nfrerun/stats")
merged_folder <- paste0(baseFolder, "20161205_Ago1KO_slam-pulse_R2/merged_1811_nfrerun/stats")
ox_folder <- paste0(baseFolder, "20160419_S2_3h-xChg_OxUnOx_pulse-chase/nextflow/1811_rerun/ox/stats")
unox_folder <- paste0(baseFolder, "20160419_S2_3h-xChg_OxUnOx_pulse-chase/nextflow/1811_rerun/unox/stats")
ago1_ip_folder <- paste0(baseFolder, "20170417_ago2ko-ago1ip-r2/merged_nf_1811_rerun/stats")
col_input_folder <- paste0(baseFolder, "20180522_S2_Ago2KO_column-input_R1_R2/1811_nextflow/01_input/stats")
col_ago_folder <- paste0(baseFolder, "20180522_S2_Ago2KO_column-input_R1_R2/1811_nextflow/02_ago-frac/stats")
col_salt_folder <- paste0(baseFolder, "20180522_S2_Ago2KO_column-input_R1_R2/1811_nextflow/03_high-salt")


filtPosFiles <- function(x) {return(paste0(x, "/filteredPositions.tsv"))}
filtLenDisFiles <- function(x) {return(paste0(x, "/filteredLenDis.tsv"))}
mutFiles <- function(x) {return(paste0(x, "/miRs.wAllMuts.tsv"))}

file_list <- c(rep1_folder, rep2_folder, merged_folder, ox_folder, unox_folder, ago1_ip_folder,
               col_input_folder, col_ago_folder)

experiment_types <- c("ago2ko-rep1", "ago2ko-rep2", "ago2ko-merge", "wt-ox", "wt-unox", "ago1-ip-merge",
                      "column-input", "column-ago-fract")

# filt.counts.File <- as_tibble(list("file" = c(paste0(rep1_folder, "/filteredPositions.tsv"),
#                                               paste0(rep2_folder, "/filteredPositions.tsv"),
#                                               paste0(merged_folder, "/filteredPositions.tsv"),
#                                               paste0(ox_folder, "/filteredPositions.tsv"),
#                                               paste0(unox_folder, "/filteredPositions.tsv"),
#                                               paste0(ago1_ip_folder, "/filteredPositions.tsv")),
#                                    "type" = c("ago2ko-rep1", "ago2ko-rep2", "ago2ko-merge",
#                                               "wt-ox", "wt-unox", "ago1-ip-merge")))
filt.counts.File <- as_tibble(list("file" = unlist(lapply(file_list, filtPosFiles)),
                                   "type" = experiment_types))
filt.lendis.File  <- as_tibble(list("file" = unlist(lapply(file_list, filtLenDisFiles)),
                                    "type" = experiment_types))
muts.File <- as_tibble(list("file" = unlist(lapply(file_list, mutFiles)),
                            "type" = experiment_types))

rep1_tp <- c(45493:45501)
rep2_tp <- c(47117:47125)
merge_tp <- c(45493:45501)
ox_tp <- c(38519:38526)
unox_tp <- c(38509:38517)
ago1_ip_tp <- c(50845:50853)
col_input_tp <- c(66842:66845, 69455, 66847:66853)
col_ago_tp <- c(67699:67710)

expDF <-
  as_tibble(list("experiment" = c(rep("ago2ko-rep1", length(rep1_tp)),
                                  rep("ago2ko-rep2", length(rep2_tp)),
                                  rep("ago2ko-merge", length(merge_tp)),
                                  rep("wt-ox", length(ox_tp)),
                                  rep("wt-unox", length(unox_tp)),
                                  rep("ago1-ip-merge", length(ago1_ip_tp)),
                                  rep("column-input", length(col_input_tp)),
                                  rep("column-ago-fract", length(col_ago_tp))),
                 "timepoint" = c(rep1_tp, rep2_tp, merge_tp, ox_tp, unox_tp, ago1_ip_tp,
                                 col_input_tp, col_ago_tp)))
readFiles <- function(file, type) { read_tsv(file) %>% mutate(experiment = type) }


filt.counts <-
  filt.counts.File %>%
  purrr::pmap(readFiles) %>%
  purrr::reduce(bind_rows) %>%
  select(-(`5p`:`3p`))
filt.lendis <-
  filt.lendis.File %>%
  purrr::pmap(readFiles) %>%
  purrr::reduce(bind_rows) %>%
  select(-(`5p`:`3p`))

# We're only filtering for miRs that we know will cause troubles down the line
# Also taking out 'wt-ox' for overlapping positions, as this experimentally is
# designed to capture different species than all others
final.filter.pos <-
  filt.counts %>%
  filter(!grepl("2a|2b-1|13b-1|276b|281-1", arm.name),
         !(experiment %in% c("wt-ox", "ago1-ip-merge")),
         !grepl("column", experiment)) %>%
  select(pos, arm.name, experiment, timepoint, time, average.ppm) %>%
  distinct() %>%
  group_by(pos, arm.name) %>%
  # Filter miR's that are present in _ALL_ libraries, except `wt-ox` and `ago1-ip`
  filter(n() == length(c(rep1_tp, rep2_tp, merge_tp, unox_tp))) %>%
  ungroup() %>%
  select(pos, arm.name) %>%
  distinct()

final.filter.counts <-
  final.filter.pos %>% 
  left_join(filt.counts) %>%
  mutate(reads = ifelse(read.type == "tcReads",
                        ifelse(UCount == 0,
                               0,
                               reads / UCount),
                        reads))
final.filter.lendis <-
  final.filter.pos %>%
  left_join(filt.lendis) %>%
  mutate(reads = ifelse(LD.type == "tcLenDis",
                        ifelse(UCount == 0,
                               0,
                               reads / UCount),
                        reads))

bg.counts <- final.filter.counts %>% filter(time == 0)

tc.bg.counts <-
  bg.counts %>%
  filter(read.type == "tcReads") %>%
  select(pos, arm.name, experiment, tc.bg.reads = reads)

bg.ratio <-
  bg.counts %>%
  filter(read.type == "totalReads") %>%
  select(pos, arm.name, experiment, stdy.state.reads = reads) %>%
  left_join(tc.bg.counts) %>%
  mutate(bg.ratio = tc.bg.reads / stdy.state.reads) %>%
  select(-stdy.state.reads, -tc.bg.reads)

# Apply cutoff to each timepoint
lendis.cutoff.raw.filtered <-
  final.filter.counts %>%
  select(pos, flybase_id, mir_name, arm.name, mir.type,
         experiment, timepoint, time, read.type, reads) %>%
  spread(read.type, reads) %>%
  left_join(bg.ratio) %>%
  mutate(cutoff = bg.ratio * totalReads) %>%
  select(-bg.ratio) %>%
  left_join(final.filter.lendis) %>%
  # Only keep reads of length isoforms > cutoff
  mutate(reads.filtered = ifelse(LD.type == "tcLenDis",
                                 ifelse(reads < cutoff | near(reads, cutoff),
                                        0,
                                        reads),
                                 reads))

pos.cutoff <-
  final.filter.lendis %>%
  filter(time == 0) %>%
  spread(LD.type, reads) %>%
  group_by(arm.name, pos, experiment) %>%
  mutate(full.bg = sum(tcLenDis) / sum(totalLenDis),
         local.bg = tcLenDis / totalLenDis) %>%
  ungroup() %>%
  replace_na(list(local.bg = 0)) %>%
  select(flybase_id, arm.name, pos, experiment, seqLen, full.bg, local.bg)

bg.seqLen <-
  final.filter.lendis %>%
  filter(time == 0) %>%
  spread(LD.type, reads) %>%
  select(flybase_id, arm.name, pos, seqLen, experiment, bg.tcReads = tcLenDis)


# `lendis.cutoff.raw.filtered` contains all miRs.
# Mutations are only available for miRs >= 5ppm.
# We will filter `lendis.cutoff.raw.filtered` later on.

## Mutations
muts <-
  muts.File %>%
  purrr::pmap(readFiles) %>%
  purrr::reduce(bind_rows) %>%
  select(-(`5p`:`3p`), -align, -full.seq)

muts.wFracts <-
  muts %>%
  # Filter duplicate miRs out as early as possible
  # Remove positions > 18 here too
  filter(!grepl("2a|2b-1|13b-1|276b|281-1", arm.name), relPos <= 18) %>%
  spread(nucleotide, count) %>%
  replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
  gather(nucleotide, count, A:T) %>%
  mutate(mutCode = ifelse(refNuc != nucleotide,
                          paste(refNuc, nucleotide, sep = '>'),
                          refNuc)) %>%
  group_by(flybase_id, timepoint, pos, start.pos, experiment) %>%
  mutate(depth = sum(count),
         mutFract = count / depth) %>%
  ungroup() %>%
  filter(grepl('>', mutCode)) %>%
  left_join(expDF) %>%
  select(-refNuc, -nucleotide, -count, -mir_name, -read.type)

bg.muts <-
  muts.wFracts %>%
  filter(time == 0) %>%
  select(pos, relPos, flybase_id, start.pos, experiment, mutCode, bg.mut = mutFract)

muts.noBG <-
  muts.wFracts %>%
  left_join(bg.muts) %>%
  mutate(bg.minus.mut = ifelse(mutFract >= bg.mut, mutFract - bg.mut, 0)) %>%
  select(-mutFract, -bg.mut) %>%
  filter(time > 0)

tc.only <-
  muts.noBG %>%
  filter(grepl("T>", mutCode)) %>%
  group_by(flybase_id, time, start.pos, mutCode, experiment) %>%
  mutate(relMutCount = sprintf("%02d", as.integer(rank(relPos))),
         relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  select(arm.name, mir.type, mir.type, start.pos, seed, UCount, experiment,
         timepoint, time, average.ppm, bg.minus.mut, relMut) %>%
  spread(relMut, bg.minus.mut) %>%
  arrange(mir.type, desc(average.ppm))

mir.meta <-
  final.filter.pos %>%
  left_join(final.filter.counts %>%
              select(flybase_id, pos, arm.name, mir.type,
                     UCount, seed, experiment, average.ppm) %>%
              distinct())

# Attention: Mutations are filtered >= 5ppm in pipeline
# Time to filter out the low coverage miRs from `lendis.cutoff.raw.filtered`
lendis.cutoff.filtered <-
  mir.meta %>%
  filter(average.ppm >= 5) %>%
  select(pos, arm.name) %>%
  distinct() %>%
  left_join(lendis.cutoff.raw.filtered)

fits <-
  muts.noBG %>%
  filter(mutCode == "T>C") %>%
  # mutate(time = time / 60) %>%
  group_by(flybase_id, arm.name, mir.type, start.pos, experiment) %>%
  do(fit.one = tryCatch(nlsLM(bg.minus.mut ~ Plat * (1 - exp(-k * time)),
                              start = list(Plat = 0.05, k = 0.5),
                              upper = c(Inf, Inf), lower = c(0, 0),
                              control = nls.lm.control(maxiter = 1000),
                              na.action = na.omit, data = .),
                        error = function(err) {
                          return(NULL)
                        }),
     fit.two = tryCatch(wrapnlsr(bg.minus.mut ~ (Plat * pFast * 0.01) *
                                   (1 - exp(-kFast * time)) +
                                   (Plat * (100 - pFast) * 0.01) *
                                   (1 - exp(-kSlow * time)),
                                 start = list(Plat = 0.05, kFast = 0.01,
                                              kSlow = 0.001, pFast = 50),
                                 # Single number is applied to all parameters
                                 lower = 0, upper = c(Inf, Inf, Inf, 100),
                                 trace = FALSE,
                                 data = .),
                        error = function(err) {
                          return(NULL)
                        }))

fit.info <-
  fits %>% mutate(fit.select = ifelse(is.null(fit.one) & is.null(fit.two),
                                      "none",
                                      ifelse(is.null(fit.one),
                                             "two.phase",
                                             ifelse(is.null(fit.two),
                                                    "one.phase",
                                                    "unknown")))) %>%
  select(flybase_id, arm.name, mir.type, start.pos, experiment, fit.select) %>%
  filter(fit.select != "none")

bothFits <- fits %>% filter(!is.null(fit.one), !is.null(fit.two))

mirs.wExpFits <-
  # Compare single and double exp fit to determine which is better by ANOVA
  bind_cols(bothFits, bothFits %>% do(aov = anova(.$fit.one, .$fit.two))) %>%
  tidy(aov) %>%
  ungroup() %>%
  filter(!is.na(p.value)) %>%
  # Cutoff for "better" is 0.05, Null hypothesis is that one.phase fits
  mutate(fit.select = ifelse(p.value >= 0.05, 'one.phase', 'two.phase')) %>%
  select(flybase_id, arm.name, mir.type, start.pos, experiment, fit.select) %>%
  # Join in those miRs that only have a single fit
  bind_rows(fit.info %>% filter(fit.select %in% c("one.phase", "two.phase"))) %>%
  # Join in fit information to fits dataframe
  left_join(fits) %>%
  # Join this in to meta data
  left_join(mir.meta)

# Getting top miRstars by half life
top_mirStar <-
  mirs.wExpFits %>%
  filter(mir.type == "star",
         fit.select == "one.phase",
         average.ppm >= 100,
         !(experiment %in% c("wt-ox", "ago1-ip-merge"))) %>%
  rowwise() %>%
  tidy(fit.one) %>%
  ungroup()

# Ago1-IP and wt-ox need to be treated sepeartely for top plateau as they preferentially retain Ago1 and Ago2 loaded species respectively.
# This results in there not being many "short lived" species in these libraries.
# Moreover, in wt-ox the assignment of mature/star will be different to wt-unox as this is dynamic per dataset.
# Therefore, we will simply select all `one.phase` species in both experiments and then ask for the fastest 10, using this as normalisation factor.
top_stable_miRs <-
  mirs.wExpFits %>%
  filter(fit.select == "one.phase",
         average.ppm >= 100,
         experiment %in% c("wt-ox", "ago1-ip-merge")) %>%
  rowwise() %>%
  tidy(fit.one) %>%
  ungroup()

total_top_miRs <- bind_rows(top_mirStar, top_stable_miRs)

top.plateaus <-
  total_top_miRs %>%
  group_by(experiment) %>%
  filter(term == "k") %>%
  top_n(10, estimate) %>%
  select(flybase_id, arm.name, start.pos, experiment) %>% 
  left_join(total_top_miRs %>% filter(term == "Plat")) %>%
  summarise(mean.p = mean(estimate))

tc.filtered <-
  lendis.cutoff.filtered %>% select(-seqLen) %>%
  filter(LD.type == "tcLenDis") %>%
  group_by(pos, arm.name, timepoint, experiment) %>%
  mutate(tc.read.sum = sum(reads.filtered)) %>%
  select(-reads, -reads.filtered) %>% distinct() %>%
  arrange(experiment, mir.type, desc(average.ppm), time) %>%
  ungroup() %>%
  left_join(top.plateaus %>% rename(labelling.efficiency = mean.p)) %>%
  mutate(lab.corrected.reads = tc.read.sum / labelling.efficiency)

seqLen.bg.norm <-
  final.filter.lendis %>%
  spread(LD.type, reads) %>%
  left_join(bg.seqLen) %>%
  left_join(pos.cutoff) %>%
  arrange(experiment, mir.type, desc(average.ppm), time) %>%
  mutate(bg.minus = ifelse(tcLenDis < (full.bg * totalLenDis) |
                             near(tcLenDis, full.bg * totalLenDis),
                           0,
                           tcLenDis - (full.bg * totalLenDis)),
         duplex = str_replace(arm.name, "-3p|-5p", ""))

seqLen.bg.norm.sumppm <-
  seqLen.bg.norm %>%
  group_by(pos, flybase_id, arm.name, duplex, mir.type, UCount, seed,
           average.ppm, time, experiment, timepoint) %>%
  summarise(iso.sum = sum(bg.minus),
            total.sum = sum(totalLenDis)) %>%
  left_join(top.plateaus %>% rename(labelling.efficiency = mean.p)) %>%
  mutate(le.norm = iso.sum / labelling.efficiency) %>%
  ungroup()

k.bio <-
  tc.filtered %>% filter(between(time,5,30)) %>%
  group_by(experiment, arm.name, pos) %>%
  do(fit = lm(lab.corrected.reads ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(arm.name, pos, experiment, k.bio = estimate)

k.bio.seqlenbgnorm <-
  seqLen.bg.norm.sumppm %>%
  filter(between(time, 5, 30)) %>%
  group_by(experiment, arm.name, pos) %>%
  do(fit = lm(le.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(arm.name, pos, experiment, k.bio = estimate)

### Figure data
## Figure 1
muts.wFracts %>%
  filter((time == -1 | time == 0 | time == 1440),
         experiment == "wt-unox") %>%
  group_by(flybase_id, arm.name, start.pos, mir.type, timepoint,
           time, current.ppm = totalReads, average.ppm, mutCode) %>%
  summarise(avg.mut.frac = mean(mutFract)) %>%
  ungroup() %>%
  mutate(mutCode = str_replace(mutCode, ">", "_"),
         time = recode(time,
                       "-1" = "no4SU-noIAA",
                       "0" = "no4SU+IAA",
                       "1440" = "24h4SU"),
         time = factor(time, levels = c("no4SU-noIAA",
                                        "no4SU+IAA",
                                        "24h4SU"))) %>%
  spread(mutCode, avg.mut.frac) %>%
  arrange(time, desc(average.ppm)) %>%
  write_tsv(paste0(baseOutput, "/Fig1/raw/Fig1C_mut-fracts.tsv"))

## Figure 2
# Figure 2B
muts.wFracts.and.avg <-
  muts.wFracts %>%
  filter(average.ppm >= 100, mutCode == "T>C", experiment == "ago2ko-merge") %>%
  mutate(duplex = str_replace(arm.name, "-3p|-5p", "")) %>%
  select(flybase_id, duplex, arm.name, start.pos, mir.type, UCount, seed, time,
         totalReads, average.ppm, relPos, mutFract) %>%
  group_by(arm.name, start.pos, time) %>%
  mutate(avg.tc.pct = mean(mutFract) * 100) %>%
  ungroup() %>%
  arrange(time, mir.type, desc(average.ppm), relPos)

duplex.avg <-
  muts.wFracts.and.avg %>%
  select(duplex, time, avg.tc.pct) %>%
  distinct() %>%
  group_by(duplex, time) %>%
  summarise(duplex.tc.pct = mean(avg.tc.pct)) %>%
  ungroup()

muts.wFracts.and.avg %>% left_join(duplex.avg) %>% write_tsv(paste0(baseOutput, "/Fig2/raw/Fig2B_TC_fract-percent.tsv"))
muts.wFracts.and.avg %>% left_join(duplex.avg) %>% select(-relPos, -mutFract) %>% distinct() %>% write_tsv(paste0(baseOutput, "/Fig2/raw/Fig2B_TC_duplex-percent.tsv")) 

# Figure 2B plot
# ecdf of 0, 5, 15, 30 min

# Figure 2C
tc.filtered %>%
  filter(average.ppm >= 100, experiment == "ago2ko-merge", time > 0) %>%
  mutate(LD.type = "tcReads") %>%
  left_join(k.bio) %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, UCount, seed, average.ppm,
         k.bio.ppm.p.min = k.bio, labelling.efficiency, LD.type, timepoint, time,
         tc.read.sum) %>%
  unite(merged.times, LD.type, timepoint, time) %>%
  spread(merged.times, tc.read.sum) %>%
  arrange(mir.type, desc(average.ppm)) %>%
  write_tsv(paste0(baseOutput, "/Fig2/raw/Fig2C_TCreads-over-time.tsv"))

seqLen.bg.norm.sumppm %>%
  filter(average.ppm >= 100, experiment == "ago2ko-merge", time > 0) %>%
  mutate(LD.type = "bg.subtr.tcReads") %>%
  left_join(k.bio.seqlenbgnorm) %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, UCount, seed, average.ppm,
         k.bio, labelling.efficiency, LD.type, timepoint, time,
         tc.read.sum = le.norm) %>%
  unite(merged.times, LD.type, timepoint, time, sep = ".") %>%
  spread(merged.times, tc.read.sum) %>%
  arrange(mir.type, desc(average.ppm)) %>%
  write_tsv(paste0(baseOutput, "/Fig2/newNorm/raw/Fig2C_TCreads-over-time.tsv"))
  
# Figure 2DEF
premiR.data <- read_tsv(paste0(baseOutput, "/Fig2/raw/allPropertiesMean.premiRNA.tf0.05.txt"))

premiR.avg.expr.frtail <-
  premiR.data %>%
  select(gname, chr, p5, isoFrac.mean.WT, align.GM.mean.WT, align.PM.mean.WT, tail.T.mean.WT) %>%
  mutate(frtail.avg.wt = align.PM.mean.WT / (align.PM.mean.WT + align.GM.mean.WT),
         fr.uridyl.avg.wt = tail.T.mean.WT / (align.PM.mean.WT + align.GM.mean.WT),
         duplex = str_replace(gname, "-RM", ""),
         duplex = ifelse(duplex == "ban", "bantam", duplex)) %>%
  group_by(gname) %>%
  top_n(1, isoFrac.mean.WT) %>%
  top_n(1, align.GM.mean.WT) %>%
  top_n(1, p5)

tc.filtered.w.kbio <-
  tc.filtered %>%
  filter(grepl("ago2ko", experiment)) %>%
  left_join(k.bio) %>%
  mutate(experiment = str_replace(experiment, "ago2ko-", ""))
  
tc.filt.kbio.merge.ppm <-
  tc.filtered.w.kbio %>%
  filter(experiment == "merge", time == 0) %>%
  select(flybase_id, arm.name, pos, merge.ppm = average.ppm) %>%
  left_join(
    tc.filtered.w.kbio %>%
      select(flybase_id, arm.name, pos, experiment, LD.type, k.bio) %>%
      distinct() %>%
      mutate(LD.type = "k.bio") %>%
      unite(bio.title, LD.type, experiment, sep = ".") %>%
      spread(bio.title, k.bio))

seqL.bgn.wkbio <-
  seqLen.bg.norm.sumppm %>%
  filter(grepl("ago2ko", experiment)) %>%
  left_join(k.bio.seqlenbgnorm) %>%
  mutate(experiment = str_replace(experiment, "ago2ko-", ""))

seqL.bgn.kbio.merge.ppm <-
  seqL.bgn.wkbio %>%
  filter(experiment == "merge", time == 0) %>%
  select(flybase_id, arm.name, pos, merge.ppm = average.ppm) %>%
  left_join(
    seqL.bgn.wkbio %>%
      select(flybase_id, arm.name, pos, experiment, k.bio) %>%
      distinct() %>%
      mutate(LD.type = "k.bio") %>%
      unite(bio.title, LD.type, experiment, sep = ".") %>%
      spread(bio.title, k.bio))

tc.filtered.w.kbio %>%
  left_join(tc.filt.kbio.merge.ppm) %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, UCount, seed, experiment,
         average.ppm, merge.ppm, k.bio.merge, k.bio.rep1, k.bio.rep2) %>%
  distinct() %>%
  mutate(duplex = str_replace(arm.name, "-3p|-5p", "")) %>%
  left_join(premiR.avg.expr.frtail) %>%
  select(flybase_id, arm.name, duplex, start.pos, mir.type, UCount, seed, experiment, average.ppm,
         merge.ppm, k.bio.merge, k.bio.rep1, k.bio.rep2, align.GM.mean.WT, align.PM.mean.WT,
         frtail.avg.wt, fr.uridyl.avg.wt) %>%
  write_tsv(paste0(baseOutput, "/Fig2/raw/Fig2DEF_premiR-ppm_fr.tailed.tsv"))

seqL.bgn.wkbio %>%
  left_join(seqL.bgn.kbio.merge.ppm) %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, UCount, seed, experiment,
         average.ppm, merge.ppm, k.bio.merge, k.bio.rep1, k.bio.rep2) %>%
  distinct() %>%
  mutate(duplex = str_replace(arm.name, "-3p|-5p", "")) %>%
  left_join(premiR.avg.expr.frtail) %>%
  select(flybase_id, arm.name, duplex, start.pos, mir.type, UCount, seed, experiment, average.ppm,
         merge.ppm, k.bio.merge, k.bio.rep1, k.bio.rep2, align.GM.mean.WT, align.PM.mean.WT,
         frtail.avg.wt, fr.uridyl.avg.wt) %>%
  arrange(mir.type, desc(merge.ppm)) %>%
  write_tsv(paste0(baseOutput, "/Fig2/newNorm/raw/Fig2DEF_premiR-ppm_fr.tailed.tsv"))

## Figure 3

mature.mir.abundance <-
  tc.filtered %>%
  filter(average.ppm >= 100, experiment == "ago2ko-merge", time == 0) %>%
  mutate(duplex = str_replace(arm.name, "-3p|-5p", "")) %>%
  select(arm.name, mir.type, duplex, average.ppm) %>% distinct() %>%
  group_by(duplex) %>% filter(n() == 2) %>% ungroup() %>%
  filter(mir.type == "mature") %>%
  select(duplex, mature.mir.ppm = average.ppm)

tc.filtered %>%
  left_join(k.bio) %>%
  filter(average.ppm >= 100, experiment == "ago2ko-merge", time > 0) %>%
  mutate(duplex = str_replace(arm.name, "-3p|-5p", ""),
         LD.type = "tcreads.LeNorm") %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, duplex, UCount, seed,
         average.ppm, k.bio, LD.type, timepoint, time, lab.corrected.reads) %>%
  unite(libname, LD.type, timepoint, time, sep = ".") %>%
  spread(libname, lab.corrected.reads) %>%
  group_by(duplex) %>% filter(n() == 2) %>% ungroup() %>%
  left_join(mature.mir.abundance) %>%
  arrange(mir.type, desc(mature.mir.ppm)) %>%
  select(flybase_id:average.ppm, mature.mir.ppm, k.bio, everything()) %>%
  write_tsv(paste0(baseOutput, "/Fig3/raw/Fig3_mir-mirstar_ss-ppm_tc-reads_kbio.tsv"))

seqLen.bg.norm.sumppm %>%
  left_join(k.bio.seqlenbgnorm) %>%
  filter(experiment == "ago2ko-merge", average.ppm >= 100) %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, duplex, UCount, seed,
         average.ppm, k.bio, timepoint, time, le.norm) %>%
  mutate(lib.type = "lab.norm.tc") %>%
  unite(l.type, lib.type, timepoint, time, sep = ".") %>%
  mutate(l.type = factor(l.type, levels = c("lab.norm.tc.45493.0", "lab.norm.tc.45494.5", "lab.norm.tc.45495.15",
                                            "lab.norm.tc.45496.30", "lab.norm.tc.45497.60", "lab.norm.tc.45498.180",
                                            "lab.norm.tc.45499.360", "lab.norm.tc.45500.720", "lab.norm.tc.45501.1440"))) %>%
  spread(l.type, le.norm) %>%
  group_by(duplex) %>% filter(n() == 2) %>% ungroup() %>%
  left_join(mature.mir.abundance) %>%
  arrange(mir.type, desc(mature.mir.ppm)) %>%
  select(flybase_id:average.ppm, mature.mir.ppm, k.bio, everything()) %>%
  write_tsv(paste0(baseOutput, "/Fig3/newNorm/raw/Fig3_mir-mirstar_ss-ppm_tc-reads_kbio.tsv"))

seqLen.bg.norm.sumppm %>%
  left_join(k.bio.seqlenbgnorm) %>%
  left_join(mature.mir.abundance) %>%
  filter(!is.na(mature.mir.ppm)) %>%
  filter(grepl("ago2ko", experiment), mature.mir.ppm >= 100) %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, duplex, UCount, seed,
         mature.mir.ppm, average.ppm, k.bio, experiment, time, le.norm) %>%
  mutate(lib.type = "lab.norm.tc") %>%
  unite(l.type, lib.type, time, sep = ".") %>%
  mutate(l.type = factor(l.type, levels = c("lab.norm.tc.0", "lab.norm.tc.5", "lab.norm.tc.15",
                                            "lab.norm.tc.30", "lab.norm.tc.60", "lab.norm.tc.180",
                                            "lab.norm.tc.360", "lab.norm.tc.720", "lab.norm.tc.1440"))) %>%
  spread(l.type, le.norm) %>%
  arrange(experiment, mir.type, desc(mature.mir.ppm)) %>%
  select(flybase_id:seed, experiment, average.ppm, mature.mir.ppm, k.bio, everything()) %>%
  write_tsv(paste0(baseOutput, "/Fig3/background_subtracted_new/raw/Fig3_mir-mirstar_ss-ppm_tc-reads_allReplicates.tsv"))

## Figure 4

lendis.w.merged.ppm <-
  lendis.cutoff.filtered %>%
  filter(grepl("ago2ko", experiment)) %>%
  left_join(lendis.cutoff.filtered %>%
              filter(experiment == "ago2ko-merge") %>%
              select(arm.name, mir.type, pos, merge.ppm = average.ppm) %>%
              distinct())

lendis.w.merged.ppm %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, UCount, seed, experiment, timepoint, time, average.ppm,
         merge.ppm, steadystate.reads = totalReads, LD.type, seqLen, reads.filtered) %>%
  left_join(tc.filtered %>% select(flybase_id, arm.name, start.pos = pos, mir.type,
                                   experiment, timepoint, time, lab.corrected.reads)) %>%
  mutate(LD.type = recode(LD.type, "totalLenDis" = "steadystate.ld", "tcLenDis" = "tc.ld")) %>%
  unite(lendis, LD.type, seqLen, sep = ".") %>%
  spread(lendis, reads.filtered) %>%
  rename(tcread.sum.LeNorm = lab.corrected.reads) %>%
  filter(merge.ppm >= 100) %>%
  arrange(experiment, mir.type, desc(merge.ppm), time) %>%
  write_tsv(paste0(baseOutput, "/Fig4/raw/Fig4BC_lendis_raw.tsv"))

seqLen.bg.norm %>%
  left_join(seqL.bgn.wkbio %>%
      filter(experiment == "merge", time == 0) %>%
      select(flybase_id, arm.name, pos, merge.ppm = average.ppm)) %>%
  left_join(seqLen.bg.norm.sumppm %>%
              select(flybase_id, arm.name, pos, experiment, timepoint, time, average.ppm,
                     labelling.efficiency, steadystate.reads = total.sum,
                     lab.corr.read.sum = le.norm)) %>%
  mutate(lab.corr.reads = bg.minus / labelling.efficiency) %>%
  filter(grepl("ago2ko", experiment), merge.ppm >= 100) %>%
  select(-tcLenDis, -bg.tcReads, -full.bg, -local.bg, -bg.minus) %>%
  gather(ld, reads, c("totalLenDis", "lab.corr.reads")) %>%
  mutate(ld = recode(ld, "totalLenDis" = "steadystate.ld", "lab.corr.reads" = "tc.ld")) %>%
  unite(lendis, ld, seqLen, sep = ".") %>%
  spread(lendis, reads) %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type, UCount, seed, experiment, timepoint, time,
         average.ppm, merge.ppm, everything()) %>%
  arrange(experiment, mir.type, desc(merge.ppm), time) %>%
  write_tsv(paste0(baseOutput, "/Fig4/background_subtracted_new/raw/Fig4BC_lendis_raw.tsv"))

lendis.w.merged.ppm %>%
  mutate(LD.type = recode(LD.type, "totalLenDis" = "sstate", "tcLenDis" = "tc")) %>%
  group_by(flybase_id, arm.name, pos, experiment, average.ppm, LD.type, timepoint, time) %>%
  summarise(w.avg.len = weighted.mean(seqLen, reads.filtered)) %>%
  ungroup() %>%
  mutate(exptype = str_sub(experiment, 8, 11)) %>%
  unite(wavg, LD.type, exptype, timepoint, time, sep = ".") %>%
  select(-experiment, -average.ppm) %>%
  spread(wavg, w.avg.len) %>%
  left_join(lendis.w.merged.ppm %>%
              select(pos, arm.name, flybase_id, experiment, average.ppm) %>%
              distinct() %>%
              spread(experiment, average.ppm)) %>%
  left_join(lendis.w.merged.ppm %>%
              filter(experiment == "ago2ko-merge") %>%
              select(pos, arm.name, flybase_id, mir.type) %>%
              distinct()) %>%
  select(flybase_id, arm.name, pos, mir.type, ago2ko.r1.ppm = `ago2ko-rep1`,
         ago2ko.r2.ppm = `ago2ko-rep2`, ago2ko.m.ppm = `ago2ko-merge`, everything()) %>%
  arrange(mir.type, desc(ago2ko.m.ppm)) %>%
  write_tsv(paste0(baseOutput, "/Fig4/raw/Fig4D_wavg-len.tsv"))

seqLen.bg.norm %>%
  filter(grepl("ago2ko", experiment)) %>%
  gather(ld, read.val, c("totalLenDis", "bg.minus")) %>%
  mutate(ld = recode(ld, "totalLenDis" = "sstate", "bg.minus" = "tc"),
         exptype = str_sub(experiment, 8, 11)) %>%
  group_by(flybase_id, arm.name, pos, experiment, exptype, average.ppm, ld, timepoint, time) %>%
  summarise(w.avg.len = weighted.mean(seqLen, read.val)) %>%
  ungroup() %>%
  unite(wavg, ld, exptype, timepoint, time, sep = ".") %>%
  select(-average.ppm) %>%
  left_join(seqLen.bg.norm.sumppm %>%
              filter(grepl("ago2ko", experiment)) %>%
              select(flybase_id, arm.name, pos, experiment, average.ppm) %>%
              distinct() %>%
              spread(experiment, average.ppm)) %>%
  left_join(seqLen.bg.norm.sumppm %>%
              filter(experiment == "ago2ko-merge") %>%
              select(pos, arm.name, flybase_id, mir.type) %>%
              distinct()) %>%
  select(-experiment) %>%
  spread(wavg, w.avg.len) %>%
  select(flybase_id, arm.name, pos, mir.type, ago2ko.r1.ppm = `ago2ko-rep1`,
         ago2ko.r2.ppm = `ago2ko-rep2`, ago2ko.m.ppm = `ago2ko-merge`, everything()) %>%
  arrange(mir.type, desc(ago2ko.m.ppm)) %>%
  write_tsv(paste0(baseOutput, "/Fig4/background_subtracted_new/raw/Fig4D_wavg-len.tsv"))

## Figure 5

mir.hl.onePhase <-
  mirs.wExpFits %>%
  filter(fit.select == "one.phase") %>%
  rowwise() %>% tidy(fit.one) %>%
  filter(term == "k") %>%
  ungroup() %>%
  select(-(UCount:average.ppm), -pos, -(std.error:p.value)) %>%
  spread(term, estimate) %>%
  mutate(hl.h = log(2) / k / 60)

mir.hl.twoPhase <-
  mirs.wExpFits %>%
  filter(fit.select == "two.phase", !is.null(fit.two), !is.na(average.ppm)) %>%
  rowwise() %>% tidy(fit.two) %>%
  filter(term == "kFast" | term == "kSlow") %>%
  ungroup() %>%
  select(-(UCount:average.ppm), -pos, -(std.error:p.value)) %>%
  spread(term, estimate) %>%
  mutate(hl.fast.h = log(2) / kFast / 60,
         hl.slow.h = log(2) / kSlow / 60)

ago2ko.tidy.hls <-
  mirs.wExpFits %>%
  filter(grepl("ago2ko", experiment)) %>%
  rowwise() %>% tidy(fit.one) %>% ungroup() %>%
  filter(term == "k") %>%
  mutate(hl = log(2) / estimate / 60) %>%
  select(-pos, -average.ppm, -term, -(std.error:p.value)) %>%
  rename(k = estimate)

muts.normMax <-
  muts.noBG %>%
  left_join(muts.noBG %>%
              filter(time == 1440) %>%
              select(flybase_id, arm.name, start.pos, relPos, experiment,
                     average.ppm, mutCode, max.mut = bg.minus.mut)) %>%
  filter(grepl("T>C", mutCode)) %>%
  group_by(flybase_id, arm.name, start.pos, time, experiment) %>%
  mutate(max.norm = bg.minus.mut / mean(max.mut)) %>%
  summarise(avg.tc = mean(bg.minus.mut),
            avg.maxNorm.tc = mean(max.norm)) %>%
  ungroup() %>%
  left_join(tc.filtered %>%
              select(flybase_id, arm.name, start.pos = pos, mir.type,
                     experiment, timepoint, time, totalReads,
                     lab.corrected.reads, average.ppm)) %>%
  left_join(tc.filt.kbio.merge.ppm %>% rename(start.pos = pos)) %>%
  filter(!is.na(mir.type)) %>%
  select(-lab.corrected.reads, -totalReads, -timepoint) %>%
  mutate(l.type = "tc.frac",
         time = factor(time, levels = c(0, 5, 15, 30, 60, 180, 360, 720, 1440))) %>%
  unite(lib, l.type, time, sep = "_")

ago2ko.raw.muts <-
  muts.normMax %>%
  filter(grepl("ago2ko", experiment)) %>%
  select(-avg.maxNorm.tc) %>%
  spread(lib, avg.tc)

ago2ko.maxNorm.muts <-
  muts.normMax %>%
  filter(grepl("ago2ko", experiment)) %>%
  select(-avg.tc) %>%
  spread(lib, avg.maxNorm.tc)

ago2ko.raw.muts %>%
  left_join(ago2ko.tidy.hls %>%
              select(-k, -fit.select) %>%
              mutate(experiment = str_replace(experiment, "ago2ko-", "hl.")) %>%
              spread(experiment, hl)) %>%
  select(flybase_id, arm.name, start.pos, mir.type, experiment, average.ppm, merge.ppm,
         hl.merge, hl.rep1, hl.rep2, tc.frac.5, tc.frac.15, tc.frac.30, tc.frac.60,
         tc.frac.180, tc.frac.360, tc.frac.720, tc.frac.1440) %>%
  mutate(experiment = factor(experiment, levels = c("ago2ko-rep1", "ago2ko-rep2", "ago2ko-merge"))) %>%
  arrange(experiment, mir.type, desc(merge.ppm)) %>%
  write_tsv(paste0(baseOutput, "/Fig5/raw/Fig5_ago2ko_tc.frac_hls.tsv"))

ago2ko.maxNorm.muts %>%
  left_join(ago2ko.tidy.hls %>%
              select(-k, -fit.select) %>%
              mutate(experiment = str_replace(experiment, "ago2ko-", "hl.")) %>%
              spread(experiment, hl)) %>%
  select(flybase_id, arm.name, start.pos, mir.type, experiment, average.ppm, merge.ppm,
         hl.merge, hl.rep1, hl.rep2, tc.frac.5, tc.frac.15, tc.frac.30, tc.frac.60,
         tc.frac.180, tc.frac.360, tc.frac.720, tc.frac.1440) %>%
  mutate(experiment = factor(experiment, levels = c("ago2ko-rep1", "ago2ko-rep2", "ago2ko-merge"))) %>%
  arrange(experiment, mir.type, desc(merge.ppm)) %>%
  write_tsv(paste0(baseOutput, "/Fig5/raw/Fig5_ago2ko_tc.maxNorm_hls.tsv"))

## Figure6

experiment.avg.ppm <-
  tc.filtered %>%
  select(flybase_id, arm.name, start.pos = pos, mir.type,
         experiment, average.ppm) %>%
  distinct()

ox.unox.ppm <-
  experiment.avg.ppm %>%
  filter(grepl("wt-unox", experiment)) %>%
  rename(unox.ppm = average.ppm) %>%
  left_join(experiment.avg.ppm %>%
              filter(grepl("wt-ox", experiment)) %>%
              rename(ox.ppm = average.ppm) %>%
              select(-mir.type, -experiment)) %>%
  mutate(ox.unox.ratio = ox.ppm / unox.ppm) %>%
  select(flybase_id, arm.name, start.pos, mir.type,
         ox.ppm, unox.ppm, ox.unox.ratio) %>%
  filter(unox.ppm >= 100) %>%
  arrange(ox.unox.ratio)

ox.unox.muts <-
  muts.noBG %>%
  filter(grepl("T>C", mutCode)) %>%
  group_by(flybase_id, arm.name, start.pos, time, experiment) %>%
  summarise(avg.tc = mean(bg.minus.mut)) %>%
  ungroup() %>%
  left_join(tc.filtered %>%
              select(flybase_id, arm.name, start.pos = pos, mir.type,
                     experiment, average.ppm) %>%
              distinct()) %>%
  filter(grepl("wt-", experiment)) %>%
  select(-mir.type, -average.ppm) %>%
  mutate(m.name = "tc.frac",
         experiment = str_replace(experiment, "wt-", ""))

muts.bg.maxNorm.tidy <-
  muts.noBG %>%
  filter(grepl("T>C", mutCode)) %>%
  left_join(muts.noBG %>%
              filter(time == 1440) %>%
              select(flybase_id, arm.name, start.pos, relPos, experiment,
                     average.ppm, mutCode, max.mut = bg.minus.mut)) %>%
  group_by(flybase_id, arm.name, start.pos, time, experiment) %>%
  mutate(max.norm = bg.minus.mut / mean(max.mut),
         relMutCount = sprintf("%02d", as.integer(rank(relPos))),
         relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  select(arm.name, mir.type, mir.type, start.pos, seed, UCount, experiment,
         timepoint, time, average.ppm, max.norm, relMut)

tc.only <-
  muts.noBG %>%
  filter(grepl("T>", mutCode)) %>%
  group_by(flybase_id, time, start.pos, mutCode, experiment) %>%
  mutate(relMutCount = sprintf("%02d", as.integer(rank(relPos))),
         relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  select(arm.name, mir.type, mir.type, start.pos, seed, UCount, experiment,
         timepoint, time, average.ppm, bg.minus.mut, relMut) %>%
  spread(relMut, bg.minus.mut) %>%
  arrange(mir.type, desc(average.ppm))

ox.unox.ppm %>%
  filter(!is.na(ox.unox.ratio)) %>%
  left_join(ox.unox.muts %>%
              filter(experiment == "ox") %>%
              unite(m.type, m.name, experiment, time, sep = ".") %>%
              spread(m.type, avg.tc)) %>%
  left_join(ox.unox.muts %>%
              filter(experiment == "unox") %>%
              unite(m.type, m.name, experiment, time, sep = ".") %>%
              spread(m.type, avg.tc)) %>%
  select(flybase_id, arm.name, start.pos, mir.type, ox.ppm, unox.ppm, ox.unox.ratio,
         tc.frac.unox.15, tc.frac.unox.30, tc.frac.unox.60, tc.frac.unox.180,
         tc.frac.unox.360, tc.frac.unox.720, tc.frac.unox.1440,
         tc.frac.ox.15, tc.frac.ox.30, tc.frac.ox.60, tc.frac.ox.180, tc.frac.ox.360,
         tc.frac.ox.720, tc.frac.ox.1440) %>%
  left_join(ago2ko.muts %>%
              filter(experiment == "ago2ko-merge") %>%
              select(-(mir.type:hl.rep2))) %>%
  write_tsv(paste0(baseOutput, "/Fig6/raw/Fig6_ox.unox.ratio_tc.fracs.tsv"))

muts.normMax.oxready <-
  muts.normMax %>%
  mutate(experiment = recode(experiment,
                             "ago2ko-merge" = "a2k-mrg",
                             "ago2ko-rep1" = "a2k-r1",
                             "ago2ko-rep2" = "a2k-r2",
                             "wt-unox" = "unox",
                             "wt-ox" = "ox")) %>%
  filter(grepl("a2k|ox", experiment)) %>%
  separate(lib, c("l.type", "time"),
           convert = TRUE, sep = "\\_")

ox.unox.ppm %>%
  filter(!is.na(ox.unox.ratio)) %>%
  left_join(muts.normMax.oxready %>%
              filter(experiment == "ox") %>%
              unite(lib, l.type, experiment, time, sep = ".") %>%
              select(-avg.tc, -(mir.type:k.bio.rep2)) %>%
              spread(lib, avg.maxNorm.tc)) %>%
  left_join(muts.normMax.oxready %>%
              filter(experiment == "unox") %>%
              unite(lib, l.type, experiment, time, sep = ".") %>%
              select(-avg.tc, -(mir.type:k.bio.rep2)) %>%
              spread(lib, avg.maxNorm.tc)) %>%
  left_join(muts.normMax.oxready %>%
              filter(experiment == "a2k-mrg") %>%
              unite(lib, l.type, experiment, time, sep = ".") %>%
              select(-avg.tc, -mir.type, -average.ppm, -(k.bio.merge:k.bio.rep2)) %>%
              spread(lib, avg.maxNorm.tc)) %>%
  select(flybase_id, arm.name, start.pos, mir.type, ox.ppm, unox.ppm,
         ox.unox.ratio, ago2ko.ppm = merge.ppm,
         tc.frac.unox.15, tc.frac.unox.30, tc.frac.unox.60, tc.frac.unox.180,
         tc.frac.unox.360, tc.frac.unox.720, tc.frac.unox.1440,
         tc.frac.ox.15, tc.frac.ox.30, tc.frac.ox.60, tc.frac.ox.180, tc.frac.ox.360,
         tc.frac.ox.720, tc.frac.ox.1440,
         `tc.frac.a2k-mrg.5`, `tc.frac.a2k-mrg.15`, `tc.frac.a2k-mrg.30`,
         `tc.frac.a2k-mrg.60`, `tc.frac.a2k-mrg.180`, `tc.frac.a2k-mrg.360`,
         `tc.frac.a2k-mrg.720`, `tc.frac.a2k-mrg.1440`) %>%
  write_tsv(paste0(baseOutput, "/Fig6/raw/Fig6_ox.unox.ratio_tc.maxNorm.tsv"))

ox.unox.ppm %>%
  filter(!is.na(ox.unox.ratio)) %>%
  left_join(muts.bg.maxNorm.tidy %>%
              filter(experiment == "ago2ko-merge") %>%
              select(arm.name, start.pos, a2k.ppm = average.ppm) %>%
              distinct()) %>%
  left_join(muts.bg.maxNorm.tidy %>%
              filter(experiment == "wt-ox") %>%
              mutate(relMut = str_replace(relMut, "TC_", "TCox_")) %>%
              select(arm.name, timepoint, time, relMut, max.norm) %>%
              spread(relMut, max.norm)) %>%
  write_tsv(paste0(baseOutput, "/Fig6/raw/Fig6-x1_ox.unox.ratio_ox-tc.maxNorm-timecourse.tsv"))

ox.unox.ppm %>%
  filter(!is.na(ox.unox.ratio)) %>%
  left_join(muts.bg.maxNorm.tidy %>%
              filter(experiment == "ago2ko-merge") %>%
              select(arm.name, start.pos, a2k.ppm = average.ppm) %>%
              distinct()) %>%
  left_join(muts.bg.maxNorm.tidy %>%
              filter(experiment == "wt-unox") %>%
              mutate(relMut = str_replace(relMut, "TC_", "TCunox_")) %>%
              select(arm.name, timepoint, time, relMut, max.norm) %>%
              spread(relMut, max.norm)) %>%
  write_tsv(paste0(baseOutput, "/Fig6/raw/Fig6-x2_ox.unox.ratio_unox-tc.maxNorm-timecourse.tsv"))

ox.unox.ppm %>%
  filter(!is.na(ox.unox.ratio)) %>%
  left_join(muts.bg.maxNorm.tidy %>%
              filter(experiment == "ago2ko-merge") %>%
              select(arm.name, start.pos, a2k.ppm = average.ppm) %>%
              distinct()) %>%
  left_join(muts.bg.maxNorm.tidy %>%
              filter(experiment == "ago2ko-merge") %>%
              mutate(relMut = str_replace(relMut, "TC_", "TCa2k_")) %>%
              select(arm.name, timepoint, time, relMut, max.norm) %>%
              spread(relMut, max.norm)) %>%
  write_tsv(paste0(baseOutput, "/Fig6/raw/Fig6-x3_ox.unox.ratio_a2k-tc.maxNorm-timecourse.tsv"))

## Figure S2
# S2A: mutation data exported for Fig 1C
# S2B and C: abundance 

## Figure S3

## Figure S4
# From wt-unox experiment
# Requires positional T>C-%

muts.wFracts %>%
  filter(grepl("T>C", mutCode), experiment == "wt-unox") %>%
  left_join(bg.muts) %>%
  mutate(bg.minus.mut = ifelse(mutFract >= bg.mut, mutFract - bg.mut, 0)) %>%
  select(-mutFract, -bg.mut) %>%
  group_by(flybase_id, time, start.pos, mutCode, experiment) %>%
  mutate(relMutCount = sprintf("%02d", as.integer(relPos)),
         relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", ""),
         UGroup = ifelse(UCount <= 3, "≤3U",
                         ifelse(UCount >= 8, "≥8U",
                                paste0(UCount, "U"))),
         UGroup = factor(UGroup, levels = c("≤3U", "4U", "5U", "6U", "7U", "≥8U")),
         arm.type = str_sub(arm.name, -2)) %>%
  ungroup() %>%
  select(arm.name, mir.type, mir.type, arm.type, start.pos, seed, UCount, UGroup, experiment,
         timepoint, time, current.ppm = totalReads, average.ppm, bg.minus.mut, relMut) %>%
  spread(relMut, bg.minus.mut) %>%
  arrange(mir.type, desc(average.ppm), time) %>%
  left_join(seqLen.bg.norm.sumppm %>%
              select(arm.name, mir.type, pos, experiment, timepoint, time, le.norm),
            by = c("arm.name", "mir.type", "experiment", "timepoint", "time",
                   "start.pos" = "pos")) %>%
  select(arm.name:time, tc.labNorm.reads = le.norm, everything()) %>%
  write_tsv(paste0(baseOutput, "/FigS4/raw/FigS4_unox_ppm_tc-fracs.tsv"))

## Figure S5

## Figure S6