library(tidyverse)

filt.counts.File <- list('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/filteredPositions.tsv',
                         '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M5001 M5013/raw/filteredPositions.tsv')
filt.lendis.File <- list('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/filteredLenDis.tsv',
                         '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M5001 M5013/raw/filteredLenDis.tsv')

ltc.tp <- c(45493:45501)
stc.tp <- c(56037:56047)

expDF <- as_data_frame(list("experiment" = c(rep('Ago2KO-24h', length(ltc.tp)),
                                             rep('Ago2KO-4h', length(stc.tp))),
                            'timepoint' = c(ltc.tp, stc.tp)))

readF <- function (x) { return(read_tsv(x)) }

filt.counts <- lapply(filt.counts.File, readF) %>% purrr::reduce(bind_rows) %>% left_join(expDF)
filt.lendis <- lapply(filt.lendis.File, readF) %>% purrr::reduce(bind_rows) %>% left_join(expDF)

tp.count <- length(c(ltc.tp, stc.tp))

# First filter out miRs & positions that are only represented in all (joined) libraries
# Then only keep miRs where the average.ppm across both libraries is maximal
# Finally, only keep these positions for joins
final.filter.pos <-
  filt.counts %>%
  group_by(pos, arm.name, read.type) %>%
  mutate(libcount = n()) %>%
  filter(libcount == tp.count) %>%
  ungroup() %>%
  select(pos, arm.name, average.ppm, experiment) %>%
  distinct() %>%
  group_by(pos, arm.name) %>% mutate(exp.avg = mean(average.ppm)) %>%
  group_by(arm.name) %>% filter(exp.avg == max(exp.avg)) %>% ungroup() %>%
  select(pos, arm.name) %>% distinct()

# U normalisation happens (almost) earliest possible
final.filter.counts <- final.filter.pos %>% left_join(filt.counts) %>% mutate(reads = ifelse(read.type == "tcReads", reads / UCount, reads))
final.filter.lendis <- final.filter.pos %>% left_join(filt.lendis) %>% mutate(reads = ifelse(LD.type == "tcLenDis", reads / UCount, reads))

# We split all reads at time = 0 to tcReads and totalReads, to calculate the background ratio
# This will be used as a cutoff for each length isoform
bg.counts <- final.filter.counts %>% filter(time == 0)

tc.bg.counts <-
  bg.counts %>%
  filter(read.type == "tcReads") %>%
  select(pos, arm.name, experiment, reads) %>%
  rename(tc.bg.reads = reads)

bg.ratio <-
  bg.counts %>%
  filter(read.type == "totalReads") %>%
  select(pos, arm.name, experiment, reads) %>%
  rename(stdy.state.reads = reads) %>%
  left_join(tc.bg.counts) %>%
  mutate(bg.ratio = tc.bg.reads / stdy.state.reads) %>%
  select(-stdy.state.reads, -tc.bg.reads)

# Apply cutoff to each timepoint
lenDis.w.cutoffs <-
  final.filter.counts %>%
  select(pos, flybase_id, mir_name, arm.name, read.type, mir.type, reads, experiment, timepoint) %>%
  spread(read.type, reads) %>%
  left_join(bg.ratio) %>%
  mutate(cutoff = bg.ratio * totalReads) %>%
  select(-bg.ratio) %>%
  left_join(final.filter.lendis) %>%
  select(-`5p`, -`3p`)

# Only keep reads of length isoforms > cutoff + 1e5 to avoid issues with floats
lendis.cutoff.filtered <-
  lenDis.w.cutoffs %>%
  mutate(reads.filtered = ifelse(LD.type == "totalLenDis", reads,
                                 ifelse(reads > cutoff + 1e-7, reads, 0))) # if LD.type == "tcLenDis"

# Short timecourse
short.mutsF <- '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M5001 M5013/raw/miRs.wAllMuts.tsv'
short.muts <- read_tsv(short.mutsF)

short.muts.wFracts <-
  short.muts %>%
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
  dplyr::select(-refNuc, -nucleotide, -count, -`5p`, -`3p`, -align, -full.seq, -mir_name, -read.type, -depth, -miRNAreads)

short.bg.muts <-
  short.muts.wFracts %>%
  filter(time == 0) %>%
  dplyr::rename(bg.mut = mutFract) %>%
  select(pos, relPos, flybase_id, start.pos, mutCode, bg.mut)

short.muts.noBG <-
  short.muts.wFracts %>%
  left_join(short.bg.muts) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-7, mutFract - bg.mut, 0)) %>%
  select(-mutFract, -bg.mut) %>%
  filter(time > 0)

short.tc.only <-
  short.muts.noBG %>%
  filter(grepl("T>C", mutCode)) %>%
  group_by(flybase_id, time, start.pos, mutCode) %>%
  mutate(relMutCount = sprintf("%02d", rank(relPos)),
         relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  dplyr::select(arm.name, mir.type, mir.type, start.pos, seed, UCount, timepoint, time, average.ppm, totalReads, bg.minus.mut, relMut) %>%
  spread(relMut, bg.minus.mut) %>%
  arrange(mir.type, desc(average.ppm))

short.tc.only %>% write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M5001 M5013/raw/tc.mutations.tsv')

# Long Timecourse
# First, we need to estimate labelling efficiency. To do this, we load the mutation file
mutsF <- '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/miRs.wAllMuts.tsv'
mutsSF <- '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M5001 M5013/raw/miRs.wAllMuts.tsv'
muts <- dplyr::bind_rows(read_tsv(mutsF), read_tsv(mutsSF))


muts.wFracts <-
  muts %>%
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
  left_join(expDF) %>%
  dplyr::select(-refNuc, -nucleotide, -count, -`5p`, -`3p`, -align, -full.seq, -mir_name, -read.type, -depth, -miRNAreads)

bg.muts <-
  muts.wFracts %>%
  filter(time == 0) %>%
  dplyr::rename(bg.mut = mutFract) %>%
  select(pos, relPos, flybase_id, start.pos, experiment, mutCode, bg.mut)

muts.noBG <-
  muts.wFracts %>%
  left_join(bg.muts) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-7, mutFract - bg.mut, 0)) %>%
  select(-mutFract, -bg.mut) %>%
  filter(time > 0)

tc.only <-
  muts.noBG %>%
  filter(grepl("T>", mutCode)) %>%
  group_by(flybase_id, time, start.pos, mutCode) %>%
  mutate(relMutCount = sprintf("%02d", rank(relPos)),
         relMut = paste(mutCode, relMutCount, sep = "_"),
         relMut = str_replace(relMut, ">", "")) %>%
  ungroup() %>%
  dplyr::select(arm.name, mir.type, mir.type, start.pos, seed, UCount, timepoint, time, average.ppm, bg.minus.mut, relMut) %>%
  spread(relMut, bg.minus.mut) %>%
  arrange(mir.type, desc(average.ppm))

tc.only %>% write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/tc.mutations.tsv')

muts.wFracts %>%
  filter(average.ppm >= 100, mutCode == "T>C", experiment == "Ago2KO-24h") %>%
  select(-experiment) %>%
  group_by(arm.name, start.pos, time, mutCode) %>%
  mutate(avg.mut.pct = mean(mutFract) * 100) %>%
  ungroup() %>%
  arrange(mir.type, desc(average.ppm), time) %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/tc.mutations.gt100ppm.tidy.tsv')

muts.wFracts %>%
  filter(average.ppm >= 100, mutCode == "T>C", experiment == "Ago2KO-24h") %>% select(-experiment) %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/tc.mutations.gt100ppm.tidy.tsv')

muts.wFracts.and.avg <-
  muts.wFracts %>%
  filter(average.ppm >= 100, mutCode == "T>C") %>%
  group_by(arm.name, start.pos, experiment, time, mutCode) %>%
  mutate(avg.mut.pct = mean(mutFract) * 100) %>%
  ungroup() %>%
  arrange(mir.type, desc(average.ppm), time)

foo <-
  muts.wFracts.and.avg %>%
  filter(time == 0, average.ppm >= 100, mir.type == "mature") %>%
  select(arm.name, avg.mut.pct) %>% distinct() %>% .$avg.mut.pct

bootresults <- boot(data=foo, statistic=Bmean, R=10000)
plot(bootresults)
boot.ci(bootresults, type = c('norm', 'basic', 'perc', 'bca'))

mir.meta <-
  muts.noBG %>%
  select(flybase_id, start.pos, arm.name, mir.type, average.ppm) %>%
  distinct()

### Calculating labelling efficiency
# Let's first calculate single exponential fits for all
# Formula: y = y0 + (Plat - y0) * (1 - exp(-k * x)) -- y0 = 0
# And compare with two phase fit:
# spanFast = (Plat - y0) * pFast * 0.01
# spanSlow = (Plat - y0) * (100 - pFast) * 0.01
# y = y0 + spanFast * (1 - exp(-kFast * x)) + spanSlow * (1 - exp(-kSlow * x))
# y0 = 0

fits <-
  muts.noBG %>%
  filter(mutCode == "T>C") %>%
  # mutate(time = time / 60) %>%
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

# Now we'll use mutReads to calculate the k in ppm per minute.
# We are NOT normalising to labelling efficiency, as we're allready in that ballpark.
maxmut <-
  muts.noBG %>%
  filter(mutCode == "T>C", time == 1440) %>%
  group_by(flybase_id, start.pos) %>%
  mutate(mean.mut = mean(bg.minus.mut), lab.norm.mut = top.plateau / mean.mut) %>%
  select(flybase_id:average.ppm, mean.mut, lab.norm.mut) %>%
  select(-time, -timepoint) %>%
  distinct() %>% arrange(mir.type, desc(average.ppm))
  
fits <-
  muts.noBG %>%
  filter(mutCode == "T>C", experiment == "Ago2KO-24h") %>% select(-experiment) %>%
  # left_join(maxmut) %>%
  # group_by(flybase_id, start.pos, time) %>%
  # mutate(avg.mut = mean(bg.minus.mut)) %>%
  # select(-bg.minus.mut, -relPos, -pos) %>%
  # distinct() %>% arrange(mir.type, desc(average.ppm))
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

bothFits <- fits %>% filter(!is.null(fit.two))

whichFits <-
  bind_cols(bothFits, bothFits %>% do(aov = anova(.$fit.one, .$fit.two))) %>%
  tidy(aov) %>%
  filter(!is.na(p.value)) %>%
  mutate(fit.select = ifelse(p.value >= 0.05, 'one.phase', 'two.phase')) %>%
  select(flybase_id, start.pos, fit.select) %>%
  bind_rows(fits %>%
              filter(is.null(fit.two)) %>%
              mutate(fit.select = "one.phase") %>%
              select(flybase_id, start.pos, fit.select))

selectedFits <- left_join(fits, whichFits)

mirs.wExpFits <- left_join(mir.meta, selectedFits)

mirs.wExpFits.wRates <-
  mirs.wExpFits %>%
  filter(!is.na(fit.select)) %>%
  mutate(actual.fit = ifelse(fit.select == "two.phase", fit.two, fit.one)) %>%
  rowwise() %>%
  tidy(actual.fit) %>%
  filter(term == "k" | term == "kFast") %>%
  select(-(std.error:p.value), -term) %>%
  dplyr::rename(k.decay = estimate)

top.onephase.miRstar <-
  mirs.wExpFits %>%
  filter(mir.type == "star", fit.select == "one.phase", average.ppm >= 100) %>%
  rowwise() %>%
  tidy(fit.one) %>%
  ungroup()

top.mirstar.by.hl <-
  top.onephase.miRstar %>%
  filter(term == "k", mir.type == "star", average.ppm >= 100) %>%
  top_n(10, estimate)

top.mirstar.by.plat <-
  top.onephase.miRstar %>%
  filter(term == "Plat", mir.type == "star", average.ppm >= 100) %>%
  top_n(10, estimate)

top.mirstar.by.hl %>%
  select(flybase_id, start.pos, arm.name) %>%
  left_join(muts.noBG) %>%
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

top.plateau <-
  top.mirstar.by.hl %>%
  select(flybase_id, start.pos, arm.name) %>%
  left_join(mirs.wExpFits) %>%
  rowwise() %>%
  tidy(fit.one) %>%
  filter(term == "Plat") %>%
  ungroup() %>%
  summarise(mean.p = mean(estimate)) %>% .$mean.p

# avg Plat of top10 fastest k:     -- Plat = 0.06507141
# avg Plat of top10 highest Plat:  -- Plat = 0.07019782

# Summarise filtered reads, and normalise by UCount & Labelling efficiency
tc.filtered <-
  lendis.cutoff.filtered %>% select(-seqLen) %>%
  filter(LD.type == "tcLenDis") %>%
  group_by(pos, arm.name, timepoint) %>%
  mutate(tc.read.sum = sum(reads.filtered)) %>%
  select(-reads, -reads.filtered) %>% distinct() %>%
  arrange(experiment, mir.type, desc(average.ppm), time) %>%
  ungroup() %>%
  mutate(labelling.efficiency = top.plateau,
         lab.corrected.reads = tc.read.sum / labelling.efficiency)

# Calculate biogenesis rate <= 15
k.bio <-
  tc.filtered %>%
  filter(time <= 15) %>%
  # mutate(time = time / 60) %>%
  group_by(experiment, arm.name, pos) %>%
    do(fit = lm(lab.corrected.reads ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>%
  dplyr::rename(k.bio = estimate)

k.bio.five <-
  tc.filtered %>% filter(time <= 5) %>%
  group_by(experiment, arm.name, pos) %>%
  do(fit = lm(lab.corrected.reads ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.5 = estimate)

k.bio.thirty <-
  tc.filtered %>% filter(time <= 30) %>%
  group_by(experiment, arm.name, pos) %>%
  do(fit = lm(lab.corrected.reads ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.30 = estimate)

k.bio.total <-
  left_join(k.bio.five, k.bio) %>%
  left_join(k.bio.thirty)


mut.ppm <-
  muts.wFracts %>%
  left_join(bg.muts) %>%
  filter(experiment == "Ago2KO-24h") %>% select(-experiment) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-9, mutFract - bg.mut, 0)) %>%
  select(-pos, -relPos, -mutFract, -bg.mut, -sRNAreads) %>%
  filter(mutCode == "T>C") %>%
  group_by(flybase_id, start.pos, time) %>%
  mutate(avg.mut = mean(bg.minus.mut),
         mut.ppm = avg.mut * totalReads,
         mut.ppm.lab.norm = mut.ppm / top.plateau) %>%
  select(-mutCode, -bg.minus.mut) %>%
  distinct()

k.bio.ppm <-
  mut.ppm %>%
  filter(time <= 15) %>%
  group_by(arm.name, start.pos) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>%
  dplyr::rename(k.bio.ppm = estimate)

k.bio.ppm.five <-
  mut.ppm %>% filter(time <= 5) %>%
  group_by(arm.name, start.pos) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.ppm.5 = estimate)

k.bio.ppm.thirty <-
  mut.ppm %>% filter(time <= 30) %>%
  group_by(arm.name, start.pos) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>% dplyr::rename(k.bio.ppm.30 = estimate)

k.bio.ppm.total <-
  left_join(k.bio.ppm.five, k.bio.ppm) %>%
  left_join(k.bio.ppm.thirty)

###
# Fig 2 data (ppm also for fig 3)
###

# Join in biogenesis rates, filter out genomic miR duplications, only select 24h timepoint
tc.filtered %>%
  select(-tcReads, -totalReads, -cutoff, -lab.corrected.reads, -LD.type) %>%
  filter(experiment == "Ago2KO-24h") %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  left_join(k.bio) %>%
  left_join(mirs.wExpFits.wRates) %>%
  unite(exp, experiment, timepoint, time, sep = ".") %>%
  spread(exp, tc.read.sum) %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/tc.filtered.reads.wKbio.tsv')

# Mut-based
mut.based.ppm.summary <-
  left_join(mirs.wExpFits.wRates, mut.ppm) %>%
  ungroup() %>%
  left_join(k.bio.ppm.thirty) %>%
  select(arm.name, start.pos, mir.type, seed, UCount, average.ppm, totalReads, k.bio.ppm.30,
         fit.select, k.decay, timepoint, time, avg.mut, mut.ppm, mut.ppm.lab.norm) %>%
  arrange(mir.type, desc(average.ppm))

mut.based.ppm.summary %>%
  select(-totalReads, -avg.mut, -mut.ppm) %>%
  mutate(libInfo = "tcReads") %>%
  filter(time > 0) %>%
  unite(lib, libInfo, timepoint, time, sep = ".") %>%
  spread(lib, mut.ppm.lab.norm) %>%
  write_tsv(file.path(outPath, 'merged.tcMutReads.wKbio.wide.tsv'))

mir.mirstar.bios <-
  tc.filtered %>%
  select(-tcReads, -totalReads, -cutoff, -lab.corrected.reads, -LD.type) %>%
  filter(experiment == "Ago2KO-24h") %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  left_join(k.bio) %>%
  left_join(mirs.wExpFits.wRates) %>%
  select(arm.name, mir.type, fit.select, average.ppm, k.bio, k.decay) %>%
  distinct()

mut.ppm %>%
  left_join(tc.filtered %>% select(-cutoff, -LD.type, -mir_name, -experiment),
            by = c("flybase_id", "timepoint", "time", "mir.type", "arm.name",
                   "average.ppm", "seed", "UCount", "totalReads", 'start.pos' = 'pos')) %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/frac.labelled.vs.total.reads.tsv')


bio.vs.ppm <-
  tc.filtered %>%
  filter(experiment == "Ago2KO-24h") %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  left_join(k.bio, by = c("pos", "arm.name", "experiment")) %>%
  left_join(mirs.wExpFits.wRates, by = c("flybase_id", "arm.name", "mir.type", "average.ppm", 'pos' = 'start.pos')) %>%
  left_join(k.bio.ppm, by = c('arm.name', 'pos' = 'start.pos'))

bio.vs.ppm.only <-
  bio.vs.ppm %>%
  filter(!is.na(fit.select)) %>%
  select(arm.name, mir.type, fit.select, average.ppm, k.bio, k.bio.ppm, k.decay) %>%
  distinct()

bio.vs.ppm.only %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/tc.bio.vs.fraction.bio.tsv')

bg.muts %>% left_join(final.filter.lendis %>% select(flybase_id, pos, arm.name, mir_name, mir.type) %>% distinct(), by = c('flybase_id', 'start.pos' = 'pos'))

bio.vs.ppm.total <-
  tc.filtered %>%
  filter(experiment == "Ago2KO-24h") %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  left_join(k.bio.total, by = c("pos", "arm.name", "experiment")) %>%
  left_join(mirs.wExpFits.wRates, by = c("flybase_id", "arm.name", "mir.type", "average.ppm", 'pos' = 'start.pos')) %>%
  left_join(k.bio.ppm.total, by = c('arm.name', 'pos' = 'start.pos'))

bio.vs.ppm.total.only <-
  bio.vs.ppm.total %>%
  filter(!is.na(fit.select)) %>%
  select(arm.name, mir.type, fit.select, average.ppm, k.bio.5, k.bio, k.bio.30, k.bio.ppm.5, k.bio.ppm, k.bio.ppm.30, k.decay) %>%
  distinct()

bio.vs.ppm.total.ppm.hl <-
  bio.vs.ppm.total.only %>%
  mutate(hl = log(2) / k.decay,
         average.ppm.hl = average.ppm / 2) %>%
  select(arm.name, mir.type, fit.select, average.ppm, hl, average.ppm.hl)

bio.vs.ppm.total.decay <-
  bio.vs.ppm.total.ppm.hl %>%
  gather(time, ppm, starts_with('average.ppm')) %>%
  mutate(time = ifelse(time == "average.ppm", 0, hl)) %>%
  left_join(bio.vs.ppm.total %>% filter(!is.na(fit.select)) %>%
              select(arm.name, mir.type, average.ppm) %>% distinct()) %>%
  group_by(arm.name, mir.type, average.ppm) %>%
  do(fit = lm(ppm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>%
  dplyr::rename(k.decay.to.hl = estimate)

bio.vs.ppm.fulltbl <-
  bio.vs.ppm.total.only %>%
  left_join(bio.vs.ppm.total.ppm.hl %>% select(-average.ppm.hl)) %>%
  left_join(bio.vs.ppm.total.decay) %>%
  mutate(pair = str_sub(arm.name, 1, -4)) %>%
  select(arm.name, pair, everything())

bio.vs.ppm.fulltbl %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/biogenesis.vs.decay.tsv')

bio.vs.ppm.miR.fulltbl <-
  bio.vs.ppm.fulltbl %>%
  select(-(k.bio.5:k.bio.ppm)) %>%
  rename(k.bio = k.bio.ppm.30, optimal.fit = fit.select) %>%
  filter(mir.type == "mature")

bio.vs.ppm.miRSTAR.fulltbl <-
  bio.vs.ppm.fulltbl %>%
  filter(mir.type == "star") %>%
  select(-arm.name, -mir.type, -(k.bio.5:k.bio.ppm)) %>%
  rename(miRSTAR.k.bio = k.bio.ppm.30, miRSTAR.fit = fit.select,
         miRSTAR.ppm = average.ppm, miRSTAR.k.decay = k.decay,
         miRSTAR.hl = hl, miRSTAR.k.decay.to.hl = k.decay.to.hl)

bio.vs.ppm.fulltbl.miR_miRSTAR.pairs <-
  left_join(bio.vs.ppm.miR.fulltbl, bio.vs.ppm.miRSTAR.fulltbl) %>%
  filter(average.ppm >= 100, miRSTAR.ppm >= 100)

# From clustering
fast.miRs <- c('mir-184', 'mir-275', 'mir-13b-2', 'mir-8', 'bantam', 'mir-14')
medium.miRs <- c('mir-2b-2', 'mir-282', 'mir-276a', 'mir-277', 'mir-305', 'mir-34')
slow.miRs <- c('mir-980', 'mir-965', 'mir-317', 'mir-33', 'mir-995', 'mir-306')

bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped <-
  bio.vs.ppm.fulltbl.miR_miRSTAR.pairs %>%
  mutate(class = ifelse(pair %in% fast.miRs, "fast",
                   ifelse(pair %in% medium.miRs, "medium",
                     ifelse(pair %in% slow.miRs, "slow", 'unknown'))))

totalReads <-
  final.filter.lendis %>%
  select(-`5p`, -`3p`) %>%
  filter(LD.type == "totalLenDis", experiment == "Ago2KO-24h") %>%
  group_by(experiment, arm.name, pos, time) %>%
  mutate(totalReads = sum(reads)) %>%
  ungroup() %>%
  select(-reads, -seqLen) %>%
  distinct() %>%
  arrange(mir.type, desc(average.ppm), time)

avg.muts.wFracts <-
  muts.wFracts %>%
  left_join(bg.muts) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-9, mutFract - bg.mut, 0)) %>%
  select(-pos, -relPos, -mutFract, -bg.mut, -sRNAreads) %>%
  filter(mutCode == "T>C") %>%
  group_by(flybase_id, start.pos, time) %>%
  mutate(avg.mut = mean(bg.minus.mut)) %>%
  select(flybase_id, timepoint, time, start.pos, mutCode, avg.mut, mir.type, average.ppm) %>%
  distinct() %>%
  arrange(mir.type, desc(average.ppm))

avg.muts.wFracts %>%
  left_join(totalReads) %>%
  mutate(mut.fract = avg.mut * totalReads) %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/mut.frac.vs.reads.tsv')

tc.filtered %>%
  left_join(avg.muts.wFracts, by = c("flybase_id", "mir.type", "timepoint", "average.ppm", "time", 'pos' = 'start.pos')) %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name),
         experiment == "Ago2KO-24h", mutCode == "T>C") %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/tcReads.vs.mutation.rate.tsv')

tc.filtered %>%
  left_join(avg.muts.wFracts, by = c("flybase_id", "mir.type", "timepoint", "average.ppm", "time", 'pos' = 'start.pos')) %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name),
         experiment == "Ago2KO-24h", mutCode == "T>C") %>%
  group_by(experiment, flybase_id, start.pos, time, mutCode) %>%
  mutate(relMutCount = rank(relPos, ties.method = "first"),
         relMut = sprintf("%s_%02d", mutCode, relMutCount)) %>%
  ungroup()

ago2ko.muts <-
  muts.noBG %>%
  filter(average.ppm >= 100, mutCode == "T>C") %>%
  group_by(flybase_id, start.pos, time, mutCode) %>%
  mutate(relMutCount = rank(relPos, ties.method = "first"),
         relMut = sprintf("%s_%02d", mutCode, relMutCount),
         experiment = "Ago2KO") %>%
  ungroup() %>%
  select(arm.name, start.pos, mir.type, seed, UCount, time, timepoint, totalReads, average.ppm, bg.minus.mut, relMut, experiment) %>%
  rename(average.reads = average.ppm,
         muts.bgMinus.norm = bg.minus.mut)
  
  select(-mutCode, -relMutCount, -relPos, -sRNAreads, -pos) %>%
  mutate(relMut = paste0("ago2ko.", relMut)) %>%
  spread(relMut, mutFract)

bio.vs.ppm.only %>% write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/tc.bio.vs.fraction.bio.tsv')

mirstar.bios <-
  mir.mirstar.bios %>%
  mutate(pair = str_sub(arm.name, 1, -4)) %>%
  filter(mir.type == "star") %>%
  dplyr::rename(star.ppm = average.ppm,
                star.bio = k.bio,
                star.decay = k.decay) %>%
  select(pair, fit.select, star.ppm, star.bio, star.decay)

mir.bios <-
  mir.mirstar.bios %>%
  mutate(pair = str_sub(arm.name, 1, -4)) %>%
  filter(mir.type == "mature", average.ppm >= 100) %>%
  dplyr::rename(mir.ppm = average.ppm,
                mir.bio = k.bio,
                mir.decay = k.decay) %>%
  left_join(mirstar.bios, by = 'pair') %>%
  filter(!is.na(star.decay))

# Load in expression data for Fig 2
expression.data <- read_tsv('/Volumes/ameres/Reichholf/sequencing/gro_pol2-chip_net_qs/GROseq_PolII-ChIP_NETseq_QS-comparison.tsv')

source.expression.w.miRbio <-
  expression.data %>%
  left_join(tc.filtered %>%
              filter(experiment == "Ago2KO-24h") %>%
              select(-tcReads, -totalReads, -cutoff, -lab.corrected.reads, -LD.type,
                     -labelling.efficiency, -tc.read.sum, -timepoint, -time) %>%
              distinct() %>%
              left_join(k.bio.ppm.thirty) %>%
              select(-experiment) %>%
              mutate(locus = str_replace(str_sub(arm.name, 1, -4), 'mir', 'miR'),
                     locus = ifelse(str_sub(arm.name, 1, -6) == "mir-2a", "miR-2a",
                              ifelse(str_sub(arm.name, 1, -6) == "mir-2b", "miR-2a", locus))),
            by = 'locus') %>%
  arrange(mir.type, desc(average.ppm)) %>%
  select(-matches('qs'), -matches('ucsd'), -matches('net'), -matches('sark'), -flybase_id, -mir_name) %>%
  select(locus, arm.name, pos, mir.type, UCount, seed, average.ppm, k.bio.ppm.30, everything()) %>%
  mutate(cistronic = ifelse(locus == 'miR-2a', 'miR-2a/b',
                      ifelse(locus == 'miR-275' | locus == 'miR-305', 'miR-275,miR-305',
                       ifelse(locus == 'miR-11' | locus == 'miR-998', 'miR-11,miR-998',
                        ifelse(locus == 'miR-34' | locus == 'miR-277', 'miR-34,miR-277',
                         ifelse(locus == 'miR-279' | locus == 'miR-996', 'miR-996,miR-279',
                          ifelse(locus == 'miR-12' | locus == 'miR-283' | locus == 'miR-304', 'miR-304,miR-283,miR-12',
                           ifelse(locus == 'miR-306' | locus == 'miR-9b' | locus == 'miR-9c', 'miR-306,miR-9b,miR-9c', 'no')))))))) %>%
  group_by(cistronic) %>%
  mutate(intron.core.gro.1 = ifelse(cistronic != 'no' & is.na(intron.core.gro.1), ifelse(sum(intron.core.gro.1, na.rm = TRUE) > 0, min(intron.core.gro.1, na.rm = TRUE), NA), intron.core.gro.1),
         intron.core.gro.2 = ifelse(cistronic != 'no' & is.na(intron.core.gro.2), ifelse(sum(intron.core.gro.2, na.rm = TRUE) > 0, min(intron.core.gro.2, na.rm = TRUE), NA), intron.core.gro.2),
         intron.core.pol2chip = ifelse(cistronic != 'no' & is.na(intron.core.pol2chip), ifelse(sum(intron.core.pol2chip, na.rm = TRUE) > 0, min(intron.core.pol2chip, na.rm = TRUE), NA), intron.core.pol2chip),
         intron.noProm.core.gro.1 = ifelse(cistronic != 'no' & is.na(intron.noProm.core.gro.1), ifelse(sum(intron.noProm.core.gro.1, na.rm = TRUE) > 0, min(intron.noProm.core.gro.1, na.rm = TRUE), NA), intron.noProm.core.gro.1),
         intron.noProm.core.gro.2 = ifelse(cistronic != 'no' & is.na(intron.noProm.core.gro.2), ifelse(sum(intron.noProm.core.gro.2, na.rm = TRUE) > 0, min(intron.noProm.core.gro.2, na.rm = TRUE), NA), intron.noProm.core.gro.2),
         intron.noProm.core.pol2chip = ifelse(cistronic != 'no' & is.na(intron.noProm.core.pol2chip), ifelse(sum(intron.noProm.core.pol2chip, na.rm = TRUE) > 0, min(intron.noProm.core.pol2chip, na.rm = TRUE), NA), intron.noProm.core.pol2chip)) %>%
  ungroup()

# Load in pre-miR sequencing data for Fig 2
preMirReads <- read_tsv('/Volumes/ameres/Reichholf/sequencing/mp/allPropertiesMean.premiRNA.tf0.05.txt')
preMirReadsWname <-
  preMirReads %>%
  mutate(mir.name = str_replace(gname, "-RM", ""),
         mir.name = ifelse(mir.name == "ban", "bantam", mir.name),
         mir.name = str_replace(mir.name, "mir", "miR"),
         totalReads = align.GM.mean.WT + align.PM.mean.WT,
         fr.tailed = align.PM.mean.WT / totalReads) %>%
  select(gname, mir.name, isoSeed, isoFrac.mean.WT, isoFrac.sd.WT,
         align.GM.mean.WT, align.PM.mean.WT, totalReads, fr.tailed)

source.expression.w.miRbio %>% write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/mir.source.expression.wKbio.tsv')

source.expression.w.miRbio %>%
  mutate(mir.name = str_sub(arm.name, 1, -4), mir.name = str_replace(mir.name, "mir", 'miR')) %>%
  left_join(preMirReadsWname) %>%
  write_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw/mir.source-and-preMir.expression.wKbio.tsv')

###
# Fig 3
###
tc.norm <-
  tc.filtered %>%
  filter(time == 180, mir.type == "mature") %>%
  select(mir_name, experiment, tc.read.sum) %>%
  dplyr::rename(norm = tc.read.sum)

normed.tc.reads <-
  tc.filtered %>%
  left_join(tc.norm) %>%
  mutate(norm.tc.reads = tc.read.sum / norm,
         norm.name = paste('norm.reads', timepoint, time, sep = "."),
         LD.type = "stdy.state.reads") %>%
  unite(lib.name, LD.type, timepoint, time)

# Hackey, set mir.types to the mir.types determined in 24h timecourse
norm.mir.types <-
  normed.tc.reads %>%
  filter(experiment == "Ago2KO-24h") %>%
  select(pos, arm.name, mir.type) %>% distinct()

avg.ppm.spread <-
  normed.tc.reads %>%
  mutate(experiment = str_replace(experiment, "-", ".")) %>%
  select(arm.name, pos, experiment, average.ppm) %>% distinct() %>%
  spread(experiment, average.ppm) %>%
  dplyr::rename(Ago2KO.24h.ppm = Ago2KO.24h, Ago2KO.4h.ppm = Ago2KO.4h)

k.bio.spread <-
  normed.tc.reads %>% left_join(k.bio) %>%
  select(arm.name, pos, experiment, k.bio) %>% distinct() %>%
  mutate(experiment = paste(str_replace(experiment, "-", "."), "kbio", sep = ".")) %>%
  spread(experiment, k.bio)

outPath <- '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/merged M3987 M4134/raw'
normed.tc.reads %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  select(arm.name, pos, seed, UCount, lib.name, totalReads) %>%
  left_join(norm.mir.types) %>%
  left_join(avg.ppm.spread) %>%
  spread(lib.name, totalReads) %>%
  arrange(mir.type, desc(Ago2KO.24h.ppm)) %>%
  write_tsv(file.path(outPath, 'merged.steady-state.ppm.tsv'))

normed.tc.reads %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  select(arm.name, pos, seed, UCount, norm.name, norm.tc.reads) %>%
  left_join(norm.mir.types) %>%
  left_join(k.bio.spread) %>%
  left_join(avg.ppm.spread) %>%
  spread(norm.name, norm.tc.reads) %>%
  arrange(mir.type, desc(Ago2KO.24h.ppm)) %>%
  write_tsv(file.path(outPath, 'merged.tc-ppm.norm-to-3h.tsv'))

normed.tc.reads %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  select(arm.name, pos, seed, UCount, lib.name, tc.read.sum, labelling.efficiency) %>%
  mutate(lib.name = str_replace(lib.name, "stdy.state", "tc.labelled")) %>%
  left_join(norm.mir.types) %>%
  left_join(k.bio.spread) %>%
  left_join(avg.ppm.spread) %>%
  spread(lib.name, tc.read.sum) %>%
  arrange(mir.type, desc(Ago2KO.24h.ppm)) %>%
  write_tsv(file.path(outPath, 'merged.tc-reads-over-time.tsv'))

muts.wFracts %>%
  left_join(bg.muts) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-9, mutFract - bg.mut, 0)) %>%
  select(-pos, -relPos, -mutFract, -bg.mut, -sRNAreads) %>%
  filter(mutCode == "T>C") %>%
  group_by(flybase_id, start.pos, timepoint) %>%
  mutate(avg.mut = mean(bg.minus.mut), experiment = "mut.TCreads",
         mut.ppm = avg.mut * totalReads) %>%
  ungroup() %>%
  select(-bg.minus.mut) %>% distinct() %>%
  filter(time > 0) %>%
  unite(lib, experiment, timepoint, time, sep = ".") %>%
  left_join(avg.ppm.spread) %>%
  left_join(k.bio.ppm.thirty) %>%
  rename(k.bio = k.bio.ppm.30) %>%
  select(-average.ppm, -totalReads, -mutCode, -avg.mut) %>%
  spread(lib, mut.ppm) %>%
  arrange(mir.type, desc(Ago2KO.24h.ppm)) %>%
  write_tsv(file.path(outPath, 'merged.tc-reads-over-time.long-short.wBio.tsv'))


normed.tc.reads %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  select(arm.name, pos, seed, UCount, norm.name, tc.read.sum) %>%
  left_join(norm.mir.types) %>%
  left_join(avg.ppm.spread) %>%
  spread(norm.name, tc.read.sum) %>%
  arrange(mir.type, desc(`Ago2KO-24h`)) %>%
  write_tsv(file.path(outPath, 'merged.tc-reads.ppm.tsv'))

# Weighted average length
wgt.avg.lens <-
  lendis.cutoff.filtered %>%
  filter(grepl("24h", experiment), LD.type == "tcLenDis") %>%
  mutate(tc.unorm = reads.filtered) %>%
  group_by(flybase_id, pos, arm.name, timepoint, time) %>%
    summarise(w.avg.len = weighted.mean(seqLen, reads.filtered)) %>%
  ungroup()

wgt.avg.background <-
  wgt.avg.lens %>%
  filter(!is.na(w.avg.len)) %>%
  group_by(flybase_id, pos) %>%
  filter(time == 15) %>%
  dplyr::rename(wavg.bg = w.avg.len)

wavg.delta.len <-
  wgt.avg.lens %>%
  filter(time >= 15) %>%
  left_join(wgt.avg.background) %>%
  group_by(flybase_id, pos) %>%
  mutate(wavg.bg = min(wavg.bg, na.rm = TRUE),
         delta.len = wavg.bg - w.avg.len)

nbr.pos <- c('mir-34', 'mir-317', 'mir-275', 'mir-2a-2', 'mir-2a-1', 'mir-263a',
             'mir-996', 'mir-9b', 'mir-2b-1', 'mir-2b-2', 'mir-305', 'mir-11')

nbr.neg <- c('mir-279', 'mir-1003', 'mir-283', 'mir-9c', 'mir-988', 'mir-304', 'mir-7', 'mir-184',
             'mir-306', 'mir-8', 'mir-980', 'mir-377', 'mir-14', 'mir-13b-1', 'mir-13b-2',
             'mir-33', 'mir-308', 'bantam', 'mir-252', 'mir-965', 'mir-995', 'mir-276a',
             'mir-998', 'mir-970', 'mir-282', 'mir-79')

wavg.delta.overview <-
  norm.mir.types %>%
  left_join(avg.ppm.spread) %>%
  left_join(wavg.delta.len %>% select(-wavg.bg)) %>%
  select(arm.name, pos, mir.type, Ago2KO.24h.ppm, Ago2KO.4h.ppm, timepoint, time, w.avg.len, delta.len) %>%
  mutate(duplex = str_sub(arm.name, 1, -4),
         nbr.substrate = ifelse(duplex %in% nbr.pos, 'yes', ifelse(duplex %in% nbr.neg, 'no', 'unknown'))) %>%
  select(-duplex)

wavg.len.wide <-
  wavg.delta.overview %>%
  mutate(wavg.text = "wavg.len") %>%
  unite(wavg.lib, wavg.text, timepoint, time, sep = ".") %>%
  select(-delta.len) %>%
  spread(wavg.lib, w.avg.len)

wavg.delta.wide <-
  wavg.delta.overview %>%
  mutate(wavg.text = "delta.len") %>%
  unite(wavg.lib, wavg.text, timepoint, time, sep = ".") %>%
  select(-w.avg.len) %>%
  spread(wavg.lib, delta.len)

wavg.len.wide %>%
  left_join(wavg.delta.wide) %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  write_tsv(file.path(outPath, 'merged.wavg-lens.tsv'))

### PCA of factors dictating steady state abundance
library(ggfortify)

autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped %>% select(average.ppm, k.bio, miRSTAR.k.decay.to.hl, hl)), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, colour = "class")

bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca <- bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped %>% mutate(average.ppm = log10(average.ppm), k.bio = log(k.bio), miRSTAR.k.decay.to.hl = log(miRSTAR.k.decay.to.hl), hl = log(hl)) %>% select(average.ppm, k.bio, miRSTAR.k.decay.to.hl, hl)
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, colour = "class", label = TRUE)
row.names(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca) <- bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped$pair
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, colour = "class", label = TRUE)
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, colour = "class", label = TRUE, loadings = TRUE, loadings.label = TRUE)
bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped %>% select(-optimal.fit, -miRSTAR.fit)
bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped %>% select(pair, average.ppm, k.bio, miRSTAR.k.decay.to.hl, hl, class)
pca.test <- prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca)

bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.raw <-
  as.matrix(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped %>%
              select(average.ppm, k.bio, miRSTAR.k.decay.to.hl, hl))

bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.log <-
  as.matrix(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped %>%
              mutate(average.ppm = log10(average.ppm), k.bio = log(k.bio),
                     miRSTAR.k.decay.to.hl = log(miRSTAR.k.decay.to.hl), hl = log(hl)) %>%
              select(average.ppm, k.bio, miRSTAR.k.decay.to.hl, hl))

bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.rank <-
  as.matrix(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped %>%
              mutate(average.ppm = rank(desc(average.ppm)), k.bio = rank(desc(k.bio)),
                     miRSTAR.k.decay.to.hl = rank(desc(miRSTAR.k.decay.to.hl)), hl = rank(desc(hl))) %>%
              select(average.ppm, k.bio, miRSTAR.k.decay.to.hl, hl))

row.names(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.raw) <- bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped$pair
row.names(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.log) <- bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped$pair
row.names(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.rank) <- bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped$pair

# Contribution in %
abs.load.raw <- abs(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.raw, scale = TRUE)$rotation)
abs.load.log <- abs(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.log, scale = TRUE)$rotation)
abs.load.rank <- abs(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.rank)$rotation)

sweep(abs.load.raw, 2, colSums(abs.load.raw), "/")
sweep(abs.load.log, 2, colSums(abs.load.log), "/")
sweep(abs.load.rank, 2, colSums(abs.load.rank), "/")

abs.load.raw.noppm <- abs(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.raw[,-1], scale = TRUE)$rotation)
abs.load.log.noppm <- abs(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.log[,-1], scale = TRUE)$rotation)
abs.load.rank.noppm <- abs(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.rank[,-1])$rotation)

sweep(abs.load.raw.noppm, 2, colSums(abs.load.raw.noppm), "/")
sweep(abs.load.log.noppm, 2, colSums(abs.load.log.noppm), "/")
sweep(abs.load.rank.noppm, 2, colSums(abs.load.rank.noppm), "/")

# Plots
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.raw, scale = TRUE), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, loadings = TRUE, loadings.label = TRUE, label = TRUE, colour = 'class')
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.log, scale = TRUE), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, loadings = TRUE, loadings.label = TRUE, label = TRUE, colour = 'class')
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.rank), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, loadings = TRUE, loadings.label = TRUE, label = TRUE, colour = 'class')

# No average.ppm
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.raw[,-1], scale = TRUE), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, loadings = TRUE, loadings.label = TRUE, label = TRUE, colour = 'class')
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.log[,-1], scale = TRUE), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, loadings = TRUE, loadings.label = TRUE, label = TRUE, colour = 'class')
autoplot(prcomp(bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped.pca.rank[,-1]), data = bio.vs.ppm.fulltbl.miR_miRSTAR.pairs.grouped, loadings = TRUE, loadings.label = TRUE, label = TRUE, colour = 'class')


##########################
# Recalculation of kbio using T>C mut rate * totalReads

all.muts <- dplyr::bind_rows(muts.noBG, short.muts.noBG) %>% left_join(expDF)

all.muts.tcOnly <-
  all.muts %>%
  filter(mutCode == "T>C") %>%
  mutate(mut.based.ppm = bg.minus.mut * totalReads) %>%
  group_by(experiment, time, timepoint, flybase_id, start.pos, arm.name, mir.type, average.ppm) %>% 
  mutate(avg.mut.based.ppm = mean(mut.based.ppm)) %>%
  ungroup() %>%
  arrange(experiment, mir.type, desc(average.ppm), time)

tcReadSum.and.mutPPM <-
  all.muts.tcOnly %>%
  left_join(tc.filtered %>%
              select(pos, flybase_id, mir_name, timepoint, labelling.efficiency, tc.read.sum),
            by = c("start.pos" = 'pos', "flybase_id", "timepoint"))

# k.bio.ppm.total already showed that tc-mut based biogenesis rate is fastest
mutPPM.bio <-
  tcReadSum.and.mutPPM %>%
  filter(time <= 15) %>%
  # mutate(time = time / 60) %>%
  select(experiment, arm.name, start.pos, relPos, time, avg.mut.based.ppm, mut.based.ppm, labelling.efficiency) %>%
  mutate(label.ppm = mut.based.ppm / labelling.efficiency) %>% filter(!is.na(labelling.efficiency)) %>%
  group_by(experiment, arm.name, start.pos) %>%
  do(fit = lm(label.ppm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>%
  dplyr::rename(k.bio = estimate)

normed.tc.reads %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  select(arm.name, pos, seed, UCount, lib.name, tc.read.sum, labelling.efficiency) %>%
  mutate(lib.name = str_replace(lib.name, "stdy.state", "tc.labelled")) %>%
  left_join(norm.mir.types) %>%
  left_join(k.bio.spread) %>%
  left_join(avg.ppm.spread) %>%
  spread(lib.name, tc.read.sum) %>%
  arrange(mir.type, desc(Ago2KO.24h.ppm))

tcMutPPM <-
  tcReadSum.and.mutPPM %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  select(arm.name, start.pos, mir.type, relPos, seed, UCount, average.ppm, labelling.efficiency, mutCode, experiment, timepoint, time, mut.based.ppm, avg.mut.based.ppm)

tcMutPPM %>%
  filter(experiment == "Ago2KO-24h") %>%
  mutate(avg.label.ppm = avg.mut.based.ppm / labelling.efficiency) %>%
  select(-mut.based.ppm, -avg.mut.based.ppm, -relPos, -mutCode) %>% distinct() %>%
  mutate(experiment = "tcMutPPM") %>%
  unite(lib, experiment, timepoint, time, sep = ".") %>%
  spread(lib, avg.label.ppm) %>%
  arrange(mir.type, desc(average.ppm)) %>%
  write_tsv(file.path(outPath, 'merged.tcMut-based.ppm.over.time.tsv'))

tcMutPPM %>%
  select(-mut.based.ppm, -relPos, -labelling.efficiency) %>% distinct() %>%
  mutate(mutCode = "raw_tcMutPPM") %>%
  unite(lib, mutCode, timepoint, time, sep = ".") %>%
  spread(lib, avg.mut.based.ppm) %>%
  arrange(experiment, mir.type, desc(average.ppm)) %>%
  write_tsv(file.path(outPath, 'merged.long-and-short.raw.tcMut-based.ppm.over.time.tsv'))


muts.ltc1.F <- '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3987 - Ago2KO 24h R1/raw/miRs.wAllMuts.tsv'
muts.ltc2.F <- '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M4134 - Ago2KO 24h R2/raw/miRs.wAllMuts.tsv'
muts.ltc.single <- dplyr::bind_rows(read_tsv(muts.ltc1.F), read_tsv(muts.ltc2.F))

ltc1 <- c(45493:45501)
ltc2 <- c(47117:47125)

expDF.single <- as_data_frame(list("experiment" = c(rep('Ago2KO-24h-R1', length(ltc1)),
                                                    rep('Ago2KO-24h-R2', length(ltc2))),
                                   'timepoint' = c(ltc1, ltc2)))

muts.single.wFracts <-
  muts.ltc.single %>%
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
  left_join(expDF.single) %>%
  dplyr::select(-refNuc, -nucleotide, -count, -`5p`, -`3p`, -align, -full.seq, -mir_name, -read.type, -depth, -miRNAreads)

bg.muts.single <-
  muts.single.wFracts %>%
  filter(time == 0) %>%
  dplyr::rename(bg.mut = mutFract) %>%
  select(pos, relPos, flybase_id, start.pos, experiment, mutCode, bg.mut)

muts.single.noBG <-
  muts.single.wFracts %>%
  left_join(bg.muts.single) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-7, mutFract - bg.mut, 0)) %>%
  select(-mutFract, -bg.mut) %>%
  filter(time > 0)

# Calculate labelling efficiency first:
fits.single <-
  muts.single.noBG %>%
  filter(mutCode == "T>C") %>%
  group_by(flybase_id, start.pos, experiment) %>%
  do(fit.one = tryCatch(nlsLM(bg.minus.mut ~ Plat * (1 - exp(-k * time)),
                              start = list(Plat = 1, k = 0.5),
                              upper = c(Inf, Inf), lower = c(0, 0),
                              control = nls.lm.control(maxiter = 1000),
                              na.action = na.omit, data = .),
                        error = function(err) {
                          return(NULL)
                        }),
     fit.two = tryCatch(wrapnlsr(bg.minus.mut ~ (Plat * pFast * 0.01) * (1 - exp(-kFast * time)) + (Plat * (100 - pFast) * 0.01) * (1 - exp(-kSlow * time)),
                                 start = list(Plat = 0.1, kFast = 0.65194720, kSlow = 0.43463147, pFast = 50),
                                 lower = c(0, 0, 0, 0), upper = c(Inf, Inf, Inf, 100),
                                 trace = FALSE,
                                 data = .),
                        error = function(err) {
                          return(NULL)
                        }))

bothFits.single <- fits.single %>% filter(!is.null(fit.two), !is.null(fit.one))

whichFits.single <-
  bind_cols(bothFits.single, bothFits.single %>% do(aov = anova(.$fit.one, .$fit.two))) %>%
  tidy(aov) %>%
  filter(!is.na(p.value)) %>%
  mutate(fit.select = ifelse(p.value >= 0.05, 'one.phase', 'two.phase')) %>%
  select(flybase_id, start.pos, fit.select, experiment) %>%
  bind_rows(fits.single %>%
              filter(is.null(fit.two), !is.null(fit.one)) %>%
              mutate(fit.select = "one.phase") %>%
              select(flybase_id, start.pos, fit.select, experiment))

selectedFits.single <- left_join(fits.single, whichFits.single) %>% filter(!is.na(fit.select))

mir.meta.single <-
  muts.single.noBG %>%
  select(flybase_id, start.pos, arm.name, mir.type, experiment, average.ppm) %>%
  distinct()

mirs.wExpFits.single <- left_join(mir.meta.single, selectedFits.single)

mirs.wExpFits.wRates.single <-
  mirs.wExpFits.single %>%
  filter(!is.na(fit.select)) %>%
  mutate(actual.fit = ifelse(fit.select == "two.phase", fit.two, fit.one)) %>%
  rowwise() %>%
  tidy(actual.fit) %>%
  filter(term == "k" | term == "kFast") %>%
  select(-(std.error:p.value), -term) %>%
  dplyr::rename(k.decay = estimate)

top.onephase.miRstar.single <-
  mirs.wExpFits.single %>%
  filter(mir.type == "star", fit.select == "one.phase", average.ppm >= 50) %>%
  rowwise() %>%
  tidy(fit.one) %>%
  ungroup()

top.mirstar.by.hl.single <-
  top.onephase.miRstar.single %>%
  filter(term == "k", mir.type == "star", average.ppm >= 100) %>%
  group_by(experiment) %>%
  top_n(10, estimate)

top.plateaus.single <-
  top.mirstar.by.hl.single %>%
  select(flybase_id, start.pos, arm.name, experiment) %>%
  left_join(mirs.wExpFits.single) %>%
  rowwise() %>%
  tidy(fit.one) %>%
  filter(term == "Plat") %>%
  ungroup() %>%
  group_by(experiment) %>%
  summarise(mean.p = mean(estimate))

mut.single.ppm <-
  muts.single.wFracts %>%
  left_join(bg.muts.single) %>%
  # filter(experiment == "Ago2KO-24h") %>% select(-experiment) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-9, mutFract - bg.mut, 0)) %>%
  select(-pos, -relPos, -mutFract, -bg.mut, -sRNAreads) %>%
  filter(mutCode == "T>C") %>%
  left_join(top.plateaus.single %>% rename(labelling.efficiency = mean.p)) %>%
  group_by(flybase_id, start.pos, time, experiment) %>%
  mutate(avg.mut = mean(bg.minus.mut),
         mut.ppm = avg.mut * totalReads,
         mut.ppm.lab.norm = mut.ppm / labelling.efficiency) %>%
  select(-mutCode, -bg.minus.mut) %>%
  distinct()

k.bio.ppm.single <-
  mut.single.ppm %>%
  filter(time <= 30) %>%
  group_by(arm.name, start.pos, experiment) %>%
  do(fit = lm(mut.ppm.lab.norm ~ 0 + time, data = .)) %>% tidy(fit) %>%
  ungroup() %>%
  select(-term, -(std.error:p.value)) %>%
  dplyr::rename(k.bio.ppm = estimate)

mut.single.ppm.spread <-
  mut.single.ppm %>%
  ungroup() %>%
  select(arm.name, start.pos, mir.type, experiment, average.ppm) %>%
  distinct() %>%
  mutate(experiment = paste0(str_replace_all(experiment, "-", "_"), ".ppm")) %>%
  group_by(arm.name, start.pos) %>%
  mutate(exp.ppm = mean(average.ppm)) %>%
  ungroup() %>%
  left_join(k.bio.ppm.single %>%
              mutate(experiment = paste0(str_replace_all(experiment, "-", "_"), ".kbio")) %>%
              spread(experiment, k.bio.ppm)) %>%
  spread(experiment, average.ppm) %>%
  arrange(mir.type, desc(exp.ppm))
