## Supplemental Figures: Negative and Positive Controls

library(tidyverse)
library(phyloseq)
library(microViz)
library(ggpubr)

#load("data/ps-decontam-filtered-counts.RData")

source("R/filter-contaminants.R")

## ---- positive control versus zymo ----

## use dataframe "decon"

# remove treatments and rename positive controls
df <- decon %>% 
  filter(!Treatment %in% c("Control", "T2", "T3", "T4")) %>% 
  mutate(Treatment = factor(Treatment, ordered = TRUE,
                            levels = c("PCP1", "PCP2", "PCP3", "Zymo"),
                            labels = c("PC1", "PC2", "PC3", "Known")))

# plot
ggplot(df, aes(x = Treatment, y = Abundance)) +
  geom_bar(aes(fill = Genus), stat = "identity", position = "fill") +
  labs(x = NULL, y = "Relative Abundance") +
  theme_pubr(legend = "right")

# save
ggsave("R/figures/pos-control-zymo.png", plot = last_plot(), dpi = 600)

## ---- negative controls ----


