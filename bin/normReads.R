library(Biostrings)
library(tidyverse)

pMFF <- '/Volumes/ameres/Reichholf/sequencing/hairpin.fa'
preMirFasta <- readDNAStringSet(pMFF)
preMirTbl <- as_tibble(list("flybase_id" = names(preMirFasta), "full.seq" = paste(preMirFasta)))

# Long TC
# Rep1
inFold <- '/Volumes/ameres/Reichholf/sequencing/20161102_Ago2KO_slam-pulse/align/bowtie_ext-bs-all/stats'
normTP <- 45498
# Rep2
# inFold <- '/Volumes/ameres/Reichholf/sequencing/20161205_Ago1KO_slam-pulse_R2/merge_aligned/stats'
# normTP <- 47122
# Short TC
# Rep 1
# inFold <- '/Volumes/ameres/Reichholf/sequencing/20170922_S2_Ago2KO_short-tc_noVirus/aligned/bowtie-bs-ext/stats'
# normTP <- 56045
# Rep2
# inFold <- '/Volumes/ameres/Reichholf/sequencing/20171005_Ago2KO_short-tc_noVirus_R2/aligned/bowtie-bs-ext/stats'
# normTP <- 56283

tPCF <- file.path(inFold, 'topPositionCounts.tsv')
aCF <- file.path(inFold, 'allCounts.tsv')

tPC <- read_tsv(tPCF)

aC <- read_tsv(aCF)

topPosWseed <-
  tPC %>% select(flybase_id, pos, mir.type) %>%
  left_join(preMirTbl) %>%
  mutate(seed = str_sub(full.seq, pos, pos + 7),
         mirBody = str_sub(full.seq, pos, pos + mirBodyLength - 1),
         UCount = str_count(mirBody, "T")) %>%
  select(-mirBody, -full.seq) %>% distinct()

allLD <-
  aC %>% select(-matches('Reads\\.')) %>%
  gather(type, reads, matches("LenDis")) %>%
  separate(type, c("LD.type", "timepoint"), convert = TRUE) %>%
  left_join(tPC %>% select(flybase_id, pos, timepoint, mir.type, average.reads) %>% distinct()) %>%
  filter(!is.na(mir.type))

bgLD <-
  allLD %>% filter(timepoint == min(timepoint)) %>%
  select(pos, seqLen, flybase_id, LD.type, reads) %>%
  replace_na(list(reads = 0))

allReadsBGMinus <-
  allLD %>% left_join(bgLD %>% rename(bg.reads = reads)) %>%
  mutate(bg.subtract = ifelse(reads - bg.reads > 0, reads - bg.reads, 0)) %>%
  group_by(pos, flybase_id, LD.type, timepoint) %>%
  mutate(read.sum = sum(bg.subtract, na.rm = TRUE)) %>%
  select(-seqLen, -reads, -bg.reads, -bg.subtract) %>%
  distinct() %>% ungroup()

maxReads <-
  allReadsBGMinus %>%
  group_by(pos, flybase_id, LD.type, mir.type) %>%
  filter(timepoint == normTP, mir.type == "mature") %>%
  ungroup() %>% select(flybase_id, LD.type, read.sum) %>% distinct() %>%
  rename(max.reads = read.sum)

# maxReadsWseed <-
#   maxReads %>%
#   left_join(topPosWseed) %>%
#   mutate(max.reads.Unorm = max.reads / UCount)

allReadsBGMinus %>%
  filter(LD.type != "totalLenDis") %>%
  unite(lendis, LD.type, timepoint, sep = ".") %>%
  left_join(topPosWseed) %>%
  spread(lendis, read.sum) %>%
  write_tsv(file.path(inFold, 'bgSubtractedReads.tsv'))

bgMinusNormedReads <-
  allReadsBGMinus %>%
  left_join(maxReads) %>%
  filter(LD.type != "totalLenDis") %>%
  left_join(topPosWseed) %>%
  mutate(max.reads.Unorm = max.reads / UCount,
         read.Unorm = read.sum / UCount,
         read.norm = read.Unorm / max.reads.Unorm) %>%
  select(-read.sum, -read.Unorm, -max.reads) %>%
  unite(lendis, LD.type, timepoint, sep = ".") %>%
  select(arm.name, pos, seed, UCount, mir.type, average.reads, lendis, read.norm) %>%
  filter(average.reads >= 1) %>%
  arrange(mir.type, desc(average.reads)) %>%
  spread(lendis, read.norm)

spreadReads <-
  tPC %>%
  filter(average.reads >= 1) %>%
  unite(read.time, read.type, timepoint, sep = ".") %>%
  left_join(topPosWseed) %>%
  select(arm.name, pos, seed, UCount, mir.type, average.reads, read.time, reads) %>%
  arrange(mir.type, desc(average.reads)) %>%
  spread(read.time, reads)

spreadReads %>%
  filter(mir.type == "mature") %>%
  select(-mir.type) %>%
  write_tsv(file.path(inFold, 'matureReadsPPM.tsv'))

spreadReads %>%
  filter(mir.type == "star") %>%
  select(-mir.type) %>%
  write_tsv(file.path(inFold, 'starReadsPPM.tsv'))

bgMinusTCreads <-
  allReadsBGMinus %>%
  filter(average.reads > 1) %>%
  filter(LD.type != "totalLenDis") %>%
  filter(timepoint > min(timepoint)) %>%
  unite(lendis, LD.type, timepoint, sep = ".") %>%
  left_join(topPosWseed) %>%
  select(arm.name, pos, seed, UCount, mir.type, average.reads, lendis, read.sum) %>%
  arrange(mir.type, desc(average.reads)) %>%
  spread(lendis, read.sum)
  
bgMinusTCreads %>%
  filter(mir.type == "mature") %>%
  select(-mir.type) %>%
  write_tsv(file.path(inFold, 'matureTCreadsBGminus.tsv'))

bgMinusTCreads %>%
  filter(mir.type == "star") %>%
  select(-mir.type) %>%
  write_tsv(file.path(inFold, 'starTCreadsBGminus.tsv'))

bgMinusNormedReads %>%
  filter(mir.type == "mature") %>%
  select(-mir.type) %>%
  write_tsv(file.path(inFold, 'bgMinusNormMature.tsv'))

bgMinusNormedReads %>%
  filter(mir.type == "star") %>%
  select(-mir.type) %>%
  write_tsv(file.path(inFold, 'bgMinusNormStar.tsv'))