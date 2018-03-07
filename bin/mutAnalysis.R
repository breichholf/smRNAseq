scissorLD <- read_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3987 - Ago2KO 24h R1/raw/bgMinus.raw.TcReads.tsv')

allMuts <- read_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3987 - Ago2KO 24h R1/raw/mutStats.tsv')

readMuts <- function(mf) { return(read_tsv(mf, col_types = cols(depth = col_double()))) }

getnormFactor <- function(mutlist) {
  mutsummary <-
    mutlist %>%
    filter(time == 1440, average.reads >= 100) %>%
    group_by(experiment, mutCode) %>%
    summarise(avg.mut = mean(mutFract)) %>%
    mutate(factor = 1 / avg.mut)
  
  return(mutsummary)
}

avgMutsPerTime <- function(ml) {
  avgMuts <-
    ml %>%
    group_by(flybase_id, arm.name, start.pos, mir.type, average.reads, experiment, time, timepoint, mutCode) %>%
      summarise(avg.mut = mean(mutFract)) %>%
    ungroup()
  
  return(avgMuts)
}

subtractMutBG <- function(mdf) {
  mutBG <-
    mdf %>%
    filter(time == 0) %>%
    select(flybase_id, start.pos, pos, relPos, mutCode, mutFract) %>%
    dplyr::rename(mutBg = mutFract)
  
  bgMinusMuts <-
    mdf %>%
    left_join(mutBG) %>%
    mutate(bgMinusMuts = ifelse(mutFract - mutBg > 0, mutFract - mutBg, 0))
  
  return(bgMinusMuts)
}

excludeSNPs <- function(ml, expDescription) {
  noSNP <-
    ml %>%
    group_by(flybase_id, time, pos) %>%
      filter(max(mutFract) < 0.75) %>%
    ungroup()
  
  return(noSNP %>% left_join(expDescription))
}

ml <- list('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3987 - Ago2KO 24h R1/raw/mutStats.tsv',
           '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M4134 - Ago2KO 24h R2/raw/mutStats.tsv',
           '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3282 - wildtype OXIDISED/raw/mutStats.tsv',
           '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3283 - wildtype UNOX/raw/mutStats.tsv')

am <- lapply(ml, readMuts)

noSNP <- lapply(am, excludeSNPs, expDescription = expDF)

mut.summary <- lapply(noSNP, getnormFactor) %>% purrr::reduce(bind_rows)

avg.muts <- lapply(noSNP, avgMutsPerTime) %>% purrr::reduce(bind_rows)

muts.bgMinus <- lapply(noSNP, subtractMutBG) %>% purrr::reduce(bind_rows)

ox.tp <- c(38519:38526)
unox.tp <- c(38509:38517)
ltc1.tp <- c(45493:45501)
ltc2.tp <- c(47117:47125)

expDF.sep <- as_data_frame(list("experiment" = c(rep("wt-oxidised", length(ox.tp)),
                                             rep('wt-unoxidised', length(unox.tp)),
                                             rep('Ago2KO-24h-R1', length(ltc1.tp)),
                                             rep('Ago2KO-24h-R2', length(ltc2.tp))),
                            'timepoint' = c(ox.tp, unox.tp, ltc1.tp, ltc2.tp)))

muts.bgWexp <-
  muts.bgMinus %>%
  left_join(expDF.sep)

avg.mutsWexp <-
  avg.muts %>%
  left_join(expDF.sep)

maxMuts <-
  muts.bgWexp %>%
  filter(time == 1440) %>%
  group_by(experiment, start.pos, flybase_id, arm.name, mutCode) %>%
    summarise(maxMutMedian = mean(bgMinusMuts)) %>%
  ungroup()

avg.bgMinusMuts <-
  muts.bgWexp %>%
  group_by(experiment, start.pos, flybase_id, time, mutCode) %>%
    summarise(muts.bgMinus = mean(bgMinusMuts)) %>%
  ungroup()

avg.bgMinusMuts.normed <-
  muts.bgWexp %>%
  left_join(maxMuts) %>%
  mutate(muts.bgMinus.norm = bgMinusMuts / maxMutMedian) %>%
  group_by(experiment, start.pos, flybase_id, arm.name, time, mutCode) %>%
    summarise(avg.bgMinus.muts = mean(muts.bgMinus.norm)) %>%
  ungroup()


mirmutF <- list('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3987 - Ago2KO 24h R1/raw/miRs.wAllMuts.tsv',
                '~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M4134 - Ago2KO 24h R2/raw/miRs.wAllMuts.tsv')

mirmuts <- lapply(mirmutF, read_tsv) %>% purrr::reduce(bind_rows)

mirmuts.wFrac <-
  mirmuts %>%
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
  dplyr::select(-refNuc, -nucleotide, -count, -`5p`, -`3p`, -align, -full.seq, -mir_name, -read.type, -depth, -miRNAreads) %>%
  left_join(expDF.sep)

mirbgmuts <-
  mirmuts.wFrac %>%
  filter(time == 0) %>%
  dplyr::rename(bg.mut = mutFract) %>%
  select(pos, relPos, flybase_id, start.pos, mutCode, bg.mut, experiment)

mirmuts.noBG <-
  mirmuts.wFrac %>%
  left_join(mirbgmuts) %>%
  mutate(bg.minus.mut = ifelse(mutFract > bg.mut + 1e-7, mutFract - bg.mut, 0)) %>%
  select(-mutFract, -bg.mut)

avg.mirmuts.nobg <-
  mirmuts.noBG %>%
  group_by(experiment, start.pos, flybase_id, time, mutCode) %>%
  summarise(muts.bgMinus = mean(bg.minus.mut)) %>% ungroup()

mirmuts.noBG.max <-
  mirmuts.noBG %>%
  filter(time == 1440) %>%
  group_by(experiment, start.pos, flybase_id, arm.name, mutCode) %>%
  summarise(mean.mut = mean(bg.minus.mut),
            median.mut = median(bg.minus.mut)) %>%
  ungroup()

mirmuts.noBG.normed <-
  mirmuts.noBG %>%
  left_join(mirmuts.noBG.max) %>%
  mutate(normed.mut = bg.minus.mut / mean.mut)
