library(tidyverse)
library(cowplot)

wt.muts <- read_tsv('~/Dropbox/PhD/data/sRNA SLAMseq REANALYSED/M3283 - wildtype UNOX/raw/mutStats.tsv', guess_max = 15000)

tcOnly <-
  wt.muts %>%
  filter(between(time, 0, 120), mutCode == "T>C") %>%
  select(-miRNAreads, -sRNAreads, -mutCode, -depth) %>%
  arrange(time, mir.type, desc(average.reads))

# 5ppm
ggplot(tcOnly %>% filter(average.reads >= 5) %>% mutate(time = as.factor(time)), aes(log10(mutFract), colour = time)) +
  stat_ecdf() +
  theme_cowplot() +
  labs(title = "T>C mutations, cutoff 5ppm",
       y = "cumulative rate")

# 10ppm
ggplot(tcOnly %>% filter(average.reads >= 10) %>% mutate(time = as.factor(time)), aes(log10(mutFract), colour = time)) +
  stat_ecdf() +
  theme_cowplot() +
  labs(title = "T>C mutations, cutoff 10ppm",
       y = "cumulative rate")

# 20ppm
ggplot(tcOnly %>% filter(average.reads >= 20) %>% mutate(time = as.factor(time)), aes(log10(mutFract), colour = time)) +
  stat_ecdf() +
  theme_cowplot() +
  labs(title = "T>C mutations, cutoff 20ppm",
       y = "cumulative rate")

# 25ppm
ggplot(tcOnly %>% filter(average.reads >= 25) %>% mutate(time = as.factor(time)), aes(log10(mutFract), colour = time)) +
  stat_ecdf() +
  theme_cowplot() +
  labs(title = "T>C mutations, cutoff 25ppm",
       y = "cumulative rate")

# 50ppm
ggplot(tcOnly %>% filter(average.reads >= 50) %>% mutate(time = as.factor(time)), aes(log10(mutFract), colour = time)) +
  stat_ecdf() +
  theme_cowplot() +
  labs(title = "T>C mutations, cutoff 50ppm",
       y = "cumulative rate")

# 75ppm
ggplot(tcOnly %>% filter(average.reads >= 75) %>% mutate(time = as.factor(time)), aes(log10(mutFract), colour = time)) +
  stat_ecdf() +
  theme_cowplot() +
  labs(title = "T>C mutations, cutoff 75ppm",
       y = "cumulative rate")

# 100ppm
ggplot(tcOnly %>% filter(average.reads >= 100) %>% mutate(time = as.factor(time)), aes(log10(mutFract), colour = time)) +
  stat_ecdf() +
  theme_cowplot() +
  labs(title = "T>C mutations, cutoff 100ppm",
       y = "cumulative rate")

ggplot(tcOnly %>% filter(average.reads >= 100) %>% mutate(time = as.factor(time)), aes(time, mutFract)) +
  geom_boxplot() +
  theme_cowplot() +
  labs(title = "T>C mutations, cutoff 100ppm",
       x = "Time (min)",
       y = "T>C mutations (log10)") +
  # coord_cartesian(ylim = c(0, 0.025)) +
  stat_compare_means(comparisons = list(c("0", "15"), c("0", "30"), c("0", "60")),
                     method = "kruskal.test")

tcAvgs <-
  tcOnly %>%
  group_by(flybase_id, start.pos, time) %>%
  mutate(avg.tc = mean(mutFract)) %>%
  select(-mutFract, -pos, -relPos) %>%
  distinct() %>%
  ungroup()

