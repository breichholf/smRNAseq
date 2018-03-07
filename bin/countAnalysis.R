library(Biostrings)
library(tidyverse)

acF <- list('/Volumes/ameres/Reichholf/sequencing/20161102_Ago2KO_slam-pulse/align/nextflow/results/stats/allCounts.tsv',
            '/Volumes/ameres/Reichholf/sequencing/20161205_Ago1KO_slam-pulse_R2/nextflow/results/stats/allCounts.tsv',
            '/Volumes/ameres/Reichholf/sequencing/20160419_S2_3h-xChg_OxUnOx_pulse-chase/nextflow/ox/results/stats/allCounts.tsv',
            '/Volumes/ameres/Reichholf/sequencing/20160419_S2_3h-xChg_OxUnOx_pulse-chase/nextflow/unox/results/stats/allCounts.tsv')

preMirFasta
preMirTbl

ox.tp <- c(38519:38526)
unox.tp <- c(38509:38517)
ltc1.tp <- c(45493:45501)
ltc2.tp <- c(47117:47125)

readACF <- function(fl) { return(read_tsv(fl)) }

gatherAc <- function(acL) {
  gAc <-
    acL %>%
    select(-matches("LenDis")) %>%
    gather(type, reads, matches("Reads\\.")) %>%
    separate(type, c('read.type', 'timepoint', 'time'), sep = "\\.", convert = TRUE) %>%
    replace_na(list(reads = 0))
  
  return(gAc)
}

ac.all <- lapply(acF, readACF)

gatheredAc <-
  lapply(ac.all, gatherAc) %>%
  purrr::reduce(bind_rows) %>%
  mutate(experiment = ifelse(timepoint %in% ox.tp, "wt-oxidised",
                             ifelse(timepoint %in% unox.tp, "wt-unoxidised",
                                    ifelse(timepoint %in% ltc1.tp, "Ago2KO-24h-R1",
                                           ifelse(timepoint %in% ltc2.tp, "Ago2KO-24h-R2", "other")))))


libcount <- length(c(ox.tp, unox.tp, ltc1.tp, ltc2.tp))

mirBodyLength <- 18

commonPos <-
  gatheredAc %>%
  filter(read.type == "totalReads", reads > 0) %>%
  select(-seqLen) %>% distinct() %>%
  group_by(flybase_id, pos) %>%
  mutate(group_count = n()) %>% ungroup() %>%
  filter(group_count == libcount) %>% select(-group_count)

cPos.wMeta <-
  commonPos %>%
  group_by(pos, flybase_id, experiment) %>%
  mutate(average.reads = mean(reads)) %>%
  group_by(pos, arm.name) %>% mutate(exp.avg = median(average.reads)) %>%
  group_by(arm.name) %>% mutate(arm.avg = mean(exp.avg)) %>%
  group_by(mir_name) %>%
    mutate(mir.type = ifelse(exp.avg == max(arm.avg), 'mature', 'star')) %>%
  ungroup() %>%
  left_join(preMirTbl) %>%
  mutate(seed = str_sub(full.seq, pos, pos + 7),
         mirBody = str_sub(full.seq, pos, pos + mirBodyLength - 1),
         UCount = str_count(mirBody, "T")) %>%
  select(-mirBody, -full.seq) %>% distinct() %>%
  dplyr::rename(totalReads = reads)

outPath <- "~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3987 - Ago2KO 24h R1/raw/"
cPos.wMeta %>% filter(grepl('-R1', experiment)) %>% select(-experiment, -exp.avg, -arm.avg) %>%
  write_tsv(file.path(outPath, 'commonTopPos.tsv'))

outPath <- "~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M4134 - Ago2KO 24h R2/raw/"
cPos.wMeta %>% filter(grepl('-R2', experiment)) %>% select(-experiment, -exp.avg, -arm.avg) %>%
  write_tsv(file.path(outPath, 'commonTopPos.tsv'))

outPath <- "~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3282 - wildtype OXIDISED/raw/"
cPos.wMeta %>% filter(grepl('wt-ox', experiment)) %>% select(-experiment, -exp.avg, -arm.avg) %>%
  write_tsv(file.path(outPath, 'commonTopPos.tsv'))

outPath <- "~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3283 - wildtype UNOX/raw/"
cPos.wMeta %>% filter(grepl('wt-unox', experiment)) %>% select(-experiment, -exp.avg, -arm.avg) %>%
  write_tsv(file.path(outPath, 'commonTopPos.tsv'))
