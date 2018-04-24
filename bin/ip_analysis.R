library(tidyverse)

ipTmutsF <- '/Volumes/ameres/Reichholf/sequencing/20170417_ago2ko-ago1ip-r2/merged_nextflow/results/stats/miRs.wAllMuts.tsv'
ipTcountsF <- '/Volumes/ameres/Reichholf/sequencing/20170417_ago2ko-ago1ip-r2/merged_nextflow/results/stats/miRs.wAllMuts.tsv'
ipT.filt.counts.F <- '/Volumes/ameres/Reichholf/sequencing/20170417_ago2ko-ago1ip-r2/merged_nextflow/results/stats/filteredPositions.tsv'
ipT.filt.lendis.F <- '/Volumes/ameres/Reichholf/sequencing/20170417_ago2ko-ago1ip-r2/merged_nextflow/results/stats/filteredLenDis.tsv'

ipT.tp <- c(50845:50853)

ipT.exp <- as_data_frame(list("experiment" = c(rep('Ago2KO-Ago1IP', length(ipT.tp))),
                              'timepoint' = c(ipT.tp)))

ipTmuts <- read_tsv(ipTmutsF)
ipT.filtcounts <- read_tsv(ipT.filt.counts.F) %>% mutate(experiment = "Ago2KO-Ago1IP")
ipT.filtLD <- read_tsv(ipT.filt.lendis.F) %>% mutate(experiment = "Ago2KO-Ago1IP")


# `final.filter.pos` is from post.analysis.R
ipT.filter.counts <- final.filter.pos %>% left_join(ipT.filtcounts) %>% mutate(reads = ifelse(read.type == "tcReads", reads / UCount, reads))
ipT.filter.lendis <- final.filter.pos %>% left_join(ipT.filtLD) %>% mutate(reads = ifelse(LD.type == "tcLenDis", reads / UCount, reads))

ipT.tc.filtered <-
  ipT.lendis.cutoff.filtered %>% select(-seqLen) %>%
  filter(LD.type == "tcLenDis") %>%
  group_by(pos, arm.name, timepoint) %>%
  mutate(tc.read.sum = sum(reads.filtered)) %>%
  select(-reads, -reads.filtered) %>% distinct() %>%
  arrange(experiment, mir.type, desc(average.ppm), time) %>%
  ungroup() %>%
  mutate(labelling.efficiency = top.plateau,
         lab.corrected.reads = tc.read.sum / labelling.efficiency)

ipT.muts.wFracts <-
  ipTmuts %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  spread(nucleotide, count) %>%
  replace_na(list(A = 0, C = 0, G = 0, T = 0)) %>%
  gather(nucleotide, count, A:T) %>%
  mutate(mutCode = ifelse(refNuc != nucleotide,
                          paste(refNuc, nucleotide, sep = '>'),
                          refNuc)) %>%
  group_by(flybase_id, timepoint, pos, start.pos) %>%
  mutate(depth = sum(count), mutFract = count / depth) %>%
  dplyr::filter(grepl('>', mutCode)) %>%
  ungroup() %>%
  mutate(experiment = "Ago2KO-Ago1IP") %>%
  dplyr::select(-refNuc, -nucleotide, -count, -`5p`, -`3p`, -align, -full.seq, -mir_name, -read.type, -depth, -miRNAreads)

ipT.bg.muts <-
  ipT.muts.wFracts %>%
  filter(time == 0) %>%
  dplyr::rename(bg.mut = mutFract) %>%
  select(pos, relPos, flybase_id, start.pos, experiment, mutCode, bg.mut)

ipT.muts.noBG <-
  ipT.muts.wFracts %>%
  left_join(ipT.bg.muts) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-7, mutFract - bg.mut, 0)) %>%
  select(-mutFract, -bg.mut) %>%
  filter(time > 0)

ipT.tc.only <-
  ipT.muts.noBG %>%
  filter(grepl("T>", mutCode)) %>%
  group_by(flybase_id, time, start.pos, mutCode) %>%
  mutate(relMutCount = sprintf("%02d", rank(relPos)),
         relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  dplyr::select(arm.name, mir.type, mir.type, start.pos, seed, UCount, timepoint, time, average.ppm, bg.minus.mut, relMut) %>%
  spread(relMut, bg.minus.mut) %>%
  arrange(mir.type, desc(average.ppm))

ipT.mir.meta <-
  ipT.muts.noBG %>%
  select(flybase_id, start.pos, arm.name, mir.type, average.ppm) %>%
  distinct()

ipT.fits <-
  ipT.muts.noBG %>%
  filter(mutCode == "T>C") %>% select(-experiment) %>%
  group_by(flybase_id, start.pos) %>%
  do(fit.one = nlsLM(bg.minus.mut ~ Plat * (1 - exp(-k * time)),
                     start = list(Plat = 1, k = 0.5),
                     upper = c(Inf, Inf), lower = c(0, 0),
                     control = nls.lm.control(maxiter = 1000),
                     na.action = na.omit, data = .),
     fit.two = tryCatch(wrapnlsr(bg.minus.mut ~ (Plat * pFast * 0.01) * (1 - exp(-kFast * time)) + (Plat * (100 - pFast) * 0.01) * (1 - exp(-kSlow * time)),
                                 start = list(Plat = 0.1, kFast = 0.65194720, kSlow = 0.43463147, pFast = 50),
                                 lower = c(0, 0, 0, 0), upper = c(Inf, Inf, Inf, 100),
                                 trace = FALSE,
                                 data = .),
                        error = function(err) {
                          return(NULL)
                        }))

ipT.bothFits <- ipT.fits %>% filter(!is.null(fit.two))

ipT.whichFits <-
  bind_cols(ipT.bothFits, ipT.bothFits %>% do(aov = anova(.$fit.one, .$fit.two))) %>%
  tidy(aov) %>%
  filter(!is.na(p.value)) %>%
  mutate(fit.select = ifelse(p.value >= 0.05, 'one.phase', 'two.phase')) %>%
  select(flybase_id, start.pos, fit.select) %>%
  bind_rows(ipT.fits %>%
              filter(is.null(fit.two)) %>%
              mutate(fit.select = "one.phase") %>%
              select(flybase_id, start.pos, fit.select))

ipT.selectedFits <- left_join(ipT.fits, ipT.whichFits)

ipT.mirs.wExpFits <- left_join(ipT.mir.meta, ipT.selectedFits)

ipT.mirs.wExpFits.wRates <-
  ipT.mirs.wExpFits %>%
  filter(!is.na(fit.select)) %>%
  mutate(actual.fit = ifelse(fit.select == "two.phase", fit.two, fit.one)) %>%
  rowwise() %>%
  tidy(actual.fit) %>%
  filter(term == "k" | term == "kFast") %>%
  select(-(std.error:p.value), -term) %>%
  dplyr::rename(k.decay = estimate)

ipT.top.onephase.miRstar <-
  ipT.mirs.wExpFits %>%
  filter(mir.type == "star", fit.select == "one.phase", average.ppm >= 100) %>%
  rowwise() %>%
  tidy(fit.one) %>%
  ungroup()

ipT.top.mirstar.by.hl <-
  ipT.top.onephase.miRstar %>%
  filter(term == "k", mir.type == "star", average.ppm >= 100) %>%
  top_n(10, estimate)

ipT.top.mirstar.by.plat <-
  ipT.top.onephase.miRstar %>%
  filter(term == "Plat", mir.type == "star", average.ppm >= 100) %>%
  top_n(10, estimate)

ipT.top.mirstar.by.hl %>%
  select(flybase_id, start.pos, arm.name) %>%
  left_join(ipT.muts.noBG) %>%
  filter(mutCode == "T>C") %>%
  group_by(arm.name, start.pos, time, mir.type, mutCode) %>%
  summarise(avg.mut = mean(bg.minus.mut)) %>%
  group_by(mir.type, mutCode) %>%
  do(fit = nlsLM(avg.mut ~ Plat * (1 - exp(-k * time)),
                 start = list(Plat = 1, k = 0.5),
                 upper = c(Inf, Inf), lower = c(0, 0),
                 control = nls.lm.control(maxiter = 1000),
                 na.action = na.omit, data = .)) %>%
  tidy(fit) %>%
  ungroup()

# by HL:    -- Plat = 0.06442492  --  k = 0.03775720
# by plat:  -- Plat = 0.06889482  --  k = 0.03300093

ipT.top.plateau <-
  ipT.top.mirstar.by.hl %>%
  select(flybase_id, start.pos, arm.name) %>%
  left_join(ipT.mirs.wExpFits) %>%
  rowwise() %>%
  tidy(fit.one) %>%
  filter(term == "Plat") %>%
  ungroup() %>%
  summarise(mean.p = mean(estimate)) %>% .$mean.p

ipT.mut.ppm <-
  ipT.muts.wFracts %>%
  left_join(ipT.bg.muts) %>%
  select(-experiment) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-9, mutFract - bg.mut, 0)) %>%
  select(-pos, -relPos, -mutFract, -bg.mut, -sRNAreads) %>%
  filter(mutCode == "T>C") %>%
  group_by(flybase_id, start.pos, time) %>%
  mutate(avg.mut = mean(bg.minus.mut),
         mut.ppm = avg.mut * totalReads,
         mut.ppm.lab.norm = mut.ppm / top.plateau) %>%
  select(-mutCode, -bg.minus.mut) %>%
  distinct()

ipT.k.bio.ppm.five <-
  ipT.mut.ppm %>% filter(time <= 5) %>%
  group_by(arm.name, start.pos) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.ppm.5 = estimate)

ipT.k.bio.ppm.fifteen <-
  ipT.mut.ppm %>% filter(time <= 15) %>%
  group_by(arm.name, start.pos) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.ppm.15 = estimate)

ipT.k.bio.ppm.thirty <-
  ipT.mut.ppm %>% filter(time <= 30) %>%
  group_by(arm.name, start.pos) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.ppm.30 = estimate)

ipT.k.bio.ppm.sixty <-
  ipT.mut.ppm %>% filter(time <= 60) %>%
  group_by(arm.name, start.pos) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.ppm.60 = estimate)

ipT.k.bio.ppm.180 <-
  ipT.mut.ppm %>% filter(time <= 180) %>%
  group_by(arm.name, start.pos) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.ppm.180 = estimate)

ipT.k.bio.ppm.compare <-
  left_join(ipT.k.bio.ppm.five, ipT.k.bio.ppm.fifteen) %>%
  left_join(ipT.k.bio.ppm.thirty) %>% left_join()
