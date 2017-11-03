library(Biostrings)
library(tidyverse)
library(broom)

### biogenesis kinetics - expected input: `allReadsBGMinus`
### linear fits to 0, 5, 15, 30.
linFit <- function(readDf, maxTime) {
  filterDf <-
    readDf %>%
    dplyr::filter(LD.type != "totalLenDis", average.reads >= 1, time <= !!maxTime)
  
  fitDf.miR <-
    filterDf %>%
    group_by(arm.name) %>%
    do(fit = lm(read.sum.miR ~ 0 + time, data = .)) %>%
    tidy(fit) %>%
    dplyr::select(arm.name, estimate) %>%
    rename(est.miR = estimate)
  
  fitDf.smRNA <-
    filterDf %>%
    group_by(arm.name) %>%
    do(fit = lm(read.sum.sRNA ~ 0 + time, data = .)) %>%
    tidy(fit) %>%
    dplyr::select(arm.name, estimate) %>%
    rename(est.sRNA = estimate)
  
  fitDf <- left_join(fitDf.miR, fitDf.smRNA)
  
  return(fitDf)
}

subBG <- function (allCounts, readOverview, topPosWSeed) {
  # Convert wide allCounts table to tidy format
  allLD <-
    allCounts %>%
    dplyr::select(-matches('Reads\\.')) %>%
    gather(type, reads, matches("LenDis")) %>%
    separate(type, c("LD.type", "timepoint"), convert = TRUE) %>%
    left_join(tPC %>% dplyr::select(flybase_id, pos, timepoint, mir.type, average.reads) %>% distinct()) %>%
    dplyr::filter(!is.na(mir.type)) %>% left_join(readOverview, by = c("timepoint" = "idx"))
  
  # Add miRNA-normed read col
  allLD.miRnorm <- 
    allLD %>%
    mutate(miRreads = reads * sRNAreads / miRNAreads) %>%
    rename(smRNAreads = reads)
  
  # Get Background
  bgLD <-
    allLD.miRnorm %>%
    dplyr::filter(time == 0) %>%
    dplyr::select(pos, seqLen, flybase_id, LD.type, smRNAreads, miRreads)
  
  # Subtract background for sRNAreads and miRNAreads normalised individually
  allReadsBGMinus <-
    allLD.miRnorm %>% left_join(bgLD %>% rename(bg.sRNA = smRNAreads, bg.miR = miRreads)) %>%
    replace_na(list(bg.miR = 0, bg.sRNA = 0)) %>%
    mutate(bg.sub.miR = ifelse(miRreads - bg.miR > 0, miRreads - bg.miR, 0),
           bg.sub.sRNA = ifelse(smRNAreads - bg.sRNA > 0, smRNAreads - bg.sRNA, 0)) %>%
    group_by(pos, flybase_id, LD.type, timepoint) %>%
    mutate(read.sum.miR = sum(bg.sub.miR, na.rm = TRUE),
           read.sum.sRNA = sum(bg.sub.sRNA, na.rm = TRUE)) %>%
    dplyr::select(-seqLen, -miRreads, -smRNAreads, -bg.sRNA, -bg.miR, -bg.sub.miR, -bg.sub.sRNA) %>%
    distinct() %>% ungroup()
  
  # Calculate linear fits for sRNAreads and miRNAreads from 0-15 and 0-30 min
  fit15 <- linFit(allReadsBGMinus, 15) %>% rename(kbio.miRNA.15 = est.miR, kbio.smRNA.15 = est.sRNA)
  fit30 <- linFit(allReadsBGMinus, 30) %>% rename(kbio.miRNA.30 = est.miR, kbio.smRNA.30 = est.sRNA)
  
  # Join fit15 and fit30 back in to allReadsBGMinus
  allReadsBGMinusWfits <-
    allReadsBGMinus %>%
    left_join(fit15) %>% left_join(fit30)
  
  return(allReadsBGMinusWfits)
}

mergeSeedInfo <- function(tcAllReads, topPos) {
  outTbl <-
    tcAllReads %>%
    select(arm.name, pos, flybase_id, mir.type, average.reads, matches('kbio')) %>%
    na.omit() %>% distinct() %>%
    left_join(topPos) %>%
    select(arm.name, pos, flybase_id, seed, UCount, mir.type, average.reads, matches('kbio'))
  
  return(outTbl)
}

# Long TC
# Rep1
# inFold <- '/Volumes/ameres/Reichholf/sequencing/20161102_Ago2KO_slam-pulse/align/results/stats/'
# readFile <- '/Volumes/ameres/Reichholf/sequencing/20161102_Ago2KO_slam-pulse/align/bowtie_ext-bs-all/20161102_smallRNAreads.txt'
# Rep2
# inFold <- '/Volumes/ameres/Reichholf/sequencing/20161205_Ago1KO_slam-pulse_R2/nf_results/stats/'
# readFile <- '/Volumes/ameres/Reichholf/sequencing/20161205_Ago1KO_slam-pulse_R2/merge_aligned/20161205_smallRNAreads.txt'
# Short TC
# Rep 1
inFold <- '/Volumes/ameres/Reichholf/sequencing/20170922_S2_Ago2KO_short-tc_noVirus/aligned/results/stats'
readFile <- '/Volumes/ameres/Reichholf/sequencing/20170922_S2_Ago2KO_short-tc_noVirus/aligned/bowtie-bs-ext/20170922_smallRNAreads.txt'
# Rep2
# inFold <- '/Volumes/ameres/Reichholf/sequencing/20171005_Ago2KO_short-tc_noVirus_R2/aligned/results/stats'
# readFile <- '/Volumes/ameres/Reichholf/sequencing/20171005_Ago2KO_short-tc_noVirus_R2/aligned/bowtie-bs-ext/20171005_smallRNAreads.txt'
# json <- '../json/samples.json'
# LNA transfection - 57521 = scrambled 24h labelling
# inFold <- '/Volumes/ameres/Reichholf/sequencing/vh/20171020_S2_LNA-transfection/results/stats/'
# normTP <- 57521

pMFF <- file.path(inFold, '../reference/hairpin.fa')
preMirFasta <- readDNAStringSet(pMFF)
preMirTbl <- as_tibble(list("flybase_id" = names(preMirFasta), "full.seq" = paste(preMirFasta)))

tPCF <- file.path(inFold, 'topPositionCounts.tsv')
aCF <- file.path(inFold, 'allCounts.tsv')

readOverview <- read_tsv(readFile)
tPC <- read_tsv(tPCF)
aC <- read_tsv(aCF)

tPC <- left_join(tPC, readOverview, by = c('timepoint' = 'idx'))

# Norm timepoint in minutes
# normTime <- 180

# normClass
normBy <- 'sRNAreads'
# normBy <- 'miRNAreads'

# Get seed and U-count from fasta
topPosWseed <-
  tPC %>% select(flybase_id, pos, mir.type) %>%
  left_join(preMirTbl) %>%
  mutate(seed = str_sub(full.seq, pos, pos + 7),
         mirBody = str_sub(full.seq, pos, pos + mirBodyLength - 1),
         UCount = str_count(mirBody, "T")) %>%
  select(-mirBody, -full.seq) %>% distinct()

# Convert wide allCounts table to tidy format
allLD <-
  aC %>% select(-matches('Reads\\.')) %>%
  gather(type, reads, matches("LenDis")) %>%
  separate(type, c("LD.type", "timepoint"), convert = TRUE) %>%
  left_join(tPC %>% select(flybase_id, pos, timepoint, mir.type, average.reads) %>% distinct()) %>%
  filter(!is.na(mir.type)) %>% left_join(readOverview, by = c("timepoint" = "idx"))

if (normBy == "miRNAreads") {
  allLD <- mutate(allLD, reads = reads * sRNAreads / miRNAreads)
  tPC <- mutate(tPC, reads = reads * sRNAreads / miRNAreads)
}

if (normBy == "miRNAreads") {
  bgMinusReadF <- 'bgSubtractedReads-miRNAnorm.tsv'
  maturePPMF <- 'matureReadsPPM-miRNAnorm.tsv'
  starPPMF <- 'starReadsPPM-miRNAnorm.tsv'
  matureTCF <-'matureTCreadsBGminus-miRNAnorm.tsv'
  starTCF <- 'starTCreadsBGminus-miRNAnorm.tsv'
  matureNormTCF <- 'bgMinusNormMature-miRNAnorm.tsv'
  starNormTCF <- 'bgMinusNormStar-miRNAnorm.tsv'
} else if (normBy == "sRNAreads") {
  bgMinusReadF <- 'bgSubtractedReads.tsv'
  maturePPMF <- 'matureReadsPPM.tsv'
  starPPMF <- 'starReadsPPM.tsv'
  matureTCF <-'matureTCreadsBGminus.tsv'
  starTCF <- 'starTCreadsBGminus.tsv'
  matureNormTCF <- 'bgMinusNormMature.tsv'
  starNormTCF <- 'bgMinusNormStar.tsv'
}

### Make library Steady state PPM spread for excel
spreadReads <-
  tPC %>%
  filter(average.reads >= 1) %>%
  unite(read.time, read.type, timepoint, sep = ".") %>%
  left_join(topPosWseed) %>%
  select(arm.name, pos, seed, UCount, mir.type, average.reads, read.time, reads) %>%
  arrange(mir.type, desc(average.reads)) %>%
  spread(read.time, reads)

### Write out mature and star reads to separate files for convenient copy-paste
spreadReads %>% filter(mir.type == "mature") %>% select(-mir.type) %>%
  write_tsv(file.path(inFold, maturePPMF))

spreadReads %>% filter(mir.type == "star") %>% select(-mir.type) %>%
  write_tsv(file.path(inFold, starPPMF))

allReadsBGMinus <- subBG(aC, readOverview, topPosWseed)

# Performed with settings from above.
longTC.R1.allReads <- allReadsBGMinus
longTC.R1.bio.duplex <-
  mergeSeedInfo(longTC.R1.allReads, topPosWseed) %>%
  rename(lTC.R1.avg.reads = average.reads,
         lTC.R1.kbio.miR = kbio.miRNA.15,
         lTC.R1.kbio.smRNA = kbio.smRNA.15) %>%
  select(arm.name, pos, seed, UCount, mir.type, matches('lTC'))

longTC.R2.allReads <- allReadsBGMinus
longTC.R2.bio.duplex <-
  mergeSeedInfo(longTC.R2.allReads, topPosWseed) %>%
  rename(lTC.R2.avg.reads = average.reads,
         lTC.R2.kbio.miR = kbio.miRNA.15,
         lTC.R2.kbio.smRNA = kbio.smRNA.15) %>%
  select(arm.name, pos, seed, UCount, mir.type, matches('lTC'))

shortTC.R1.allReads <- allReadsBGMinus
shortTC.R1.bio.duplex <-
  mergeSeedInfo(shortTC.R1.allReads, topPosWseed) %>%
  rename(sTC.R1.avg.reads = average.reads,
         sTC.R1.kbio.miR = kbio.miRNA.15,
         sTC.R1.kbio.smRNA = kbio.smRNA.15) %>%
  select(arm.name, pos, seed, UCount, mir.type, matches('sTC'))

shortTC.R2.allReads <- allReadsBGMinus
shortTC.R2.bio.duplex <-
  mergeSeedInfo(shortTC.R2.allReads, topPosWseed) %>%
  rename(sTC.R2.avg.reads = average.reads,
         sTC.R2.kbio.miR = kbio.miRNA.15,
         sTC.R2.kbio.smRNA = kbio.smRNA.15) %>%
  select(arm.name, pos, seed, UCount, mir.type, matches('sTC'))

### Get maximum reads AFTER background subtraction
#   normalise to `normTime` (3h)
maxReads <-
  allReadsBGMinus %>%
  group_by(pos, flybase_id, LD.type, mir.type) %>%
  filter(timepoint == normTime, mir.type == "mature") %>%
  ungroup() %>% select(flybase_id, LD.type, read.sum) %>% distinct() %>%
  rename(max.reads = read.sum)

### Write out background subtracted reads
allReadsBGMinus %>%
  filter(LD.type != "totalLenDis") %>%
  unite(lendis, LD.type, timepoint, sep = ".") %>%
  left_join(topPosWseed) %>%
  spread(lendis, read.sum) %>%
  write_tsv(file.path(inFold, bgMinusReadF))

### Make spread table for excel
bgMinusTCreads <-
  allReadsBGMinus %>%
  filter(average.reads > 1) %>%
  filter(LD.type != "totalLenDis") %>%
  filter(time > 0) %>%
  unite(lendis, LD.type, timepoint, sep = ".") %>%
  left_join(topPosWseed) %>%
  left_join(fit15) %>%
  left_join(fit30) %>%
  select(arm.name, pos, seed, UCount, mir.type, average.reads, kbio.15, kbio.30, lendis, read.sum) %>%
  arrange(mir.type, desc(average.reads)) %>%
  spread(lendis, read.sum)

### Write out mature and star reads split, and arranged descending for copy-paste
bgMinusTCreads %>% filter(mir.type == "mature") %>% select(-mir.type) %>%
  write_tsv(file.path(inFold, matureTCF))

bgMinusTCreads %>% filter(mir.type == "star") %>% select(-mir.type) %>%
  write_tsv(file.path(inFold, starTCF))

######
### Normalise background subtracted reads by U-count
######
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

### Write out mature and star reads split and arranged descinding for copy-paste
bgMinusNormedReads %>% filter(mir.type == "mature") %>% select(-mir.type) %>%
  write_tsv(file.path(inFold, matureNormTCF))

bgMinusNormedReads %>% filter(mir.type == "star") %>% select(-mir.type) %>%
  write_tsv(file.path(inFold, starNormTCF))
