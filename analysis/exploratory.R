# Data exploration
library(tidyverse)

# `tc.filtered` is a tidy tibble, with one value of tc.reads per timepoint
# `muts.noBG` are all mutations for the long timecourse
# `short.muts.noBG` -- " -- short timecourse

all.muts <- dplyr::bind_rows(muts.noBG, short.muts.noBG) %>% left_join(expDF)

all.muts.tcOnly <-
  all.muts %>%
  filter(mutCode == "T>C") %>%
  mutate(mut.based.ppm = bg.minus.mut * totalReads) %>%
  group_by(experiment, time, timepoint, flybase_id, start.pos, arm.name, mir.type, average.ppm) %>% 
  summarise(avg.mut.based.ppm = mean(mut.based.ppm)) %>%
  ungroup() %>%
  arrange(experiment, mir.type, desc(average.ppm), time)


both.tc.countedPPM_mutDerivedPPM <-
  tc.filtered %>%
  left_join(all.muts.tcOnly %>% select(-average.ppm),
            by = c("flybase_id", "arm.name", "mir.type", "experiment",
                   "timepoint", "time", 'pos' = 'start.pos')) %>%
  select(-flybase_id, -mir_name, -seed, -labelling.efficiency, -LD.type) %>%
  filter(!grepl("2a", arm.name), !grepl("2b-1", arm.name), !grepl("13b-1", arm.name), !grepl("276b", arm.name)) %>%
  replace_na(list('avg.mut.based.ppm' = 0))

both.tc.meta <-
  both.tc.countedPPM_mutDerivedPPM %>%
  filter(experiment == "Ago2KO-24h", average.ppm >= 100, time == 180) %>%
  select(pos, arm.name, mir.type)

both.tc.min.info <-
  both.tc.meta %>%
  left_join(both.tc.countedPPM_mutDerivedPPM %>% filter(time == 180) %>%
              select(pos, arm.name, mir.type, experiment, timepoint, time, average.ppm, tc.read.sum, avg.mut.based.ppm))

ggplot(both.tc.min.info %>% mutate(experiment = str_replace(experiment, "-", "_tcReadSum_")) %>%
         select(-timepoint, -time, -average.ppm, -avg.mut.based.ppm) %>% spread(experiment, tc.read.sum),
       aes(log10(Ago2KO_tcReadSum_4h), log10(Ago2KO_tcReadSum_24h))) +
  geom_point(size = 4, alpha = 0.4) +
  geom_smooth()

ggplot(both.tc.min.info %>% mutate(experiment = str_replace(experiment, "-", "_mutBased_")) %>%
         select(-timepoint, -time, -average.ppm, -tc.read.sum) %>% spread(experiment, avg.mut.based.ppm),
       aes(log10(Ago2KO_mutBased_4h), log10(Ago2KO_mutBased_24h))) + geom_point(size = 4, alpha = 0.4)

# Correlation plots

both.tc.test <-
  both.tc.countedPPM_mutDerivedPPM %>%
  select(pos, arm.name, mir.type, experiment, timepoint, time, average.ppm, totalReads, tc.read.sum, avg.mut.based.ppm) %>%
  filter(time %in% c(5, 15, 30, 60, 180)) %>%
  rename(tc.based.reads = tc.read.sum, total.based.Reads = totalReads) %>% gather(read.type, reads, matches('based'))

both.tc.test.avg.ppm <-
  both.tc.test %>%
  select(pos, arm.name, mir.type, experiment, average.ppm) %>% distinct() %>%
  mutate(experiment = str_replace(experiment, "-", "_"),
         experiment = paste0(experiment, ".avgppm")) %>%
  spread(experiment, average.ppm)

both.tc.test.long.wide <-
  both.tc.test.avg.ppm %>%
  mutate(Ago2KO_24h.avgppm = log10(Ago2KO_24h.avgppm),
         Ago2KO_4h.avgppm = log10(Ago2KO_4h.avgppm)) %>%
  left_join(both.tc.test %>%
              filter(experiment == "Ago2KO-24h" & read.type != "total.based.Reads") %>%
              select(-average.ppm, -experiment) %>%
              mutate(read.type = str_replace(read.type, '.based', ''),
                     read.type = str_replace(read.type, 'avg.mut', 'mut')) %>%
              unite(lib, read.type, timepoint, time, sep = "_") %>%
              mutate(reads = log10(reads)) %>%
              spread(lib, reads)) %>%
  arrange(mir.type, desc(Ago2KO_24h.avgppm))

both.tc.test.short.wide <-
  both.tc.test.avg.ppm %>%
  mutate(Ago2KO_24h.avgppm = log10(Ago2KO_24h.avgppm),
         Ago2KO_4h.avgppm = log10(Ago2KO_4h.avgppm)) %>%
  left_join(both.tc.test %>%
              filter(experiment == "Ago2KO-4h" & read.type != "total.based.Reads") %>%
              select(-average.ppm, -experiment) %>%
              mutate(read.type = str_replace(read.type, '.based', ''),
                     read.type = str_replace(read.type, 'avg.mut', 'mut')) %>%
              unite(lib, read.type, timepoint, time, sep = "_") %>%
              mutate(reads = log10(reads)) %>%
              spread(lib, reads)) %>%
  arrange(mir.type, desc(Ago2KO_24h.avgppm))

both.tc.test.tcRComp.wide <-
  both.tc.test.avg.ppm %>%
  mutate(Ago2KO_24h.avgppm = log10(Ago2KO_24h.avgppm),
         Ago2KO_4h.avgppm = log10(Ago2KO_4h.avgppm)) %>%
  left_join(both.tc.test %>%
              filter(read.type == "tc.based.reads") %>%
              select(-average.ppm, -read.type) %>%
              mutate(experiment = str_replace(experiment, "Ago2KO", "A2K")) %>%
              unite(lib, experiment, timepoint, time, sep = "_") %>%
              mutate(reads = log10(reads)) %>%
              spread(lib, reads)) %>%
  arrange(mir.type, desc(Ago2KO_24h.avgppm))

both.tc.test.tcMComp.wide <-
  both.tc.test.avg.ppm %>%
  mutate(Ago2KO_24h.avgppm = log10(Ago2KO_24h.avgppm),
         Ago2KO_4h.avgppm = log10(Ago2KO_4h.avgppm)) %>%
  left_join(both.tc.test %>%
              filter(read.type == "avg.mut.based.ppm") %>%
              select(-average.ppm, -read.type) %>%
              mutate(experiment = str_replace(experiment, "Ago2KO", "A2K")) %>%
              unite(lib, experiment, timepoint, time, sep = "_") %>%
              mutate(reads = log10(reads)) %>%
              spread(lib, reads)) %>%
  arrange(mir.type, desc(Ago2KO_24h.avgppm))

points_with_cor <- function(data, mapping, ..., method = "pearson") {
  x <- eval(mapping$x, data)
  y <- eval(mapping$y, data)
  cor <- cor(x, y, method = method)
  ggally_points(data, mapping, ...) +
    ggplot2::geom_label(
      data = data.frame(
        x = min(x, na.rm = TRUE),
        y = max(y, na.rm = TRUE),
        lab = round(cor, digits = 3)
      ),
      mapping = ggplot2::aes(x = x, y = y, label = lab),
      hjust = 0, vjust = 1,
      size = 5, fontface = "bold",
      inherit.aes = FALSE # do not inherit anything from the ...
    ) +
    scale_color_brewer(palette = "Set1")
}

brewer_densityDiag <- function(data, mapping, ...) {
  ggally_densityDiag(data, mapping, ...) +
    scale_fill_brewer(palette = "Set1")
}

pl <- ggpairs(both.tc.test.long.wide %>% filter(Ago2KO_24h.avgppm >= 2),
              mapping = aes(color = mir.type),
              columns = 4:ncol(both.tc.test.long.wide),
              upper = NULL,
              diag = list(continuous = wrap(brewer_densityDiag, alpha = 0.75)),
              lower = list(continuous = wrap(points_with_cor, method = "spearman", alpha = 0.45)),
              title = "Long TC only", xlab = "log10(ppm)", ylab = "log10(ppm)")

ps <- ggpairs(both.tc.test.short.wide %>% filter(Ago2KO_24h.avgppm >= 2),
              mapping = aes(color = mir.type),
              columns = 4:ncol(both.tc.test.short.wide),
              upper = NULL,
              diag = list(continuous = wrap(brewer_densityDiag, alpha = 0.75)),
              lower = list(continuous = wrap(points_with_cor, method = "spearman", alpha = 0.45)),
              title = "Short TC only", xlab = "log10(ppm)", ylab = "log10(ppm)")

pReadComp <- ggpairs(both.tc.test.tcRComp.wide %>% filter(Ago2KO_24h.avgppm >= 2),
                     mapping = aes(color = mir.type),
                     columns = 4:ncol(both.tc.test.tcRComp.wide),
                     upper = NULL,
                     diag = list(continuous = wrap(brewer_densityDiag, alpha = 0.75)),
                     lower = list(continuous = wrap(points_with_cor, method = "spearman", alpha = 0.45)),
                     title = "TC reads only", xlab = "log10(ppm)", ylab = "log10(ppm)")

pMutComp <- ggpairs(both.tc.test.tcMComp.wide %>% filter(Ago2KO_24h.avgppm >= 2),
                    mapping = aes(color = mir.type),
                    columns = 4:ncol(both.tc.test.tcMComp.wide),
                    upper = NULL,
                    diag = list(continuous = wrap(brewer_densityDiag, alpha = 0.75)),
                    lower = list(continuous = wrap(points_with_cor, method = "spearman", alpha = 0.45)),
                    title = "mut derived TC reads only", xlab = "log10(ppm)", ylab = "log10(ppm)")

ggsave('both.tc.long.pdf', pl, width = 520, height = 490, units = 'mm')
ggsave('both.tc.short.pdf', ps, width = 520, height = 490, units = 'mm')
ggsave('both.tc.tcRead_based.pdf', pReadComp, width = 520, height = 490, units = 'mm')
ggsave('both.tc.tcMut_based.pdf', pMutComp, width = 520, height = 490, units = 'mm')

