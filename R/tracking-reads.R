
## quality check tracking reads through dada2

# Emily Van Syoc 11/2021

require(tidyverse)
require(dada2)
require(rstatix)
require(ggpubr)

# get data
# all Seqkit stat output
tf <- read.table("data/reads-tracking/trimstatf.txt", sep = "\t", header = TRUE) %>% mutate(step = "trim")
#tr <- read.table("dada2-data/trimstatr.txt", sep = "\t", header = TRUE) %>% mutate(direction = "reverse")
raw <- read.table("data/reads-tracking/rawstat.txt", sep = "\t", header = TRUE) %>% mutate(step = "raw")
ad <- read.table("data/reads-tracking/adapter-remstat.txt", sep = "\t", header = TRUE) %>% mutate(step = "adapt")

# dada output
trk <- read.table("data/reads-tracking/metformintrack-reads-dada.txt", sep = "\t", header = TRUE)

## ---- wrangle ----

## wrangle

# raw files
rawf <- raw %>% 
  #mutate(direction = if_else(str_detect(file, "_1.fq.gz"), "forward", "reverse")) %>% 
  # remove reverse reads 
  filter(str_detect(file, "_1.fq.gz")) %>% 
  mutate(ID = str_extract(file, "M(\\d){1,3}")) %>% 
  select(num_seqs, ID, step)

# adapter removed files
adf <- ad %>% 
  # remove accidental extra "CM" files
  filter(!str_detect(file, "CM")) %>% 
  #mutate(direction = if_else(str_detect(file, "trimmed_1.fastq"), "forward", "reverse")) %>% 
  # remove reverse files
  filter(str_detect(file, "trimmed_1.fastq")) %>% 
  mutate(ID = str_extract(file, "M(\\d){1,3}")) %>% 
  select(num_seqs, ID, step)

# trimmed files
tff <- tf %>% 
  # get ids
  mutate(ID = str_extract(file, "M(\\d){1,3}")) %>% 
  select(num_seqs, ID, step)

# combine all seqkits
both <- rbind(rawf, adf, tff)


# make track-reads vertical
trkv <- trk  %>% 
  rownames_to_column(var = "ID") %>% 
  # remove reverse reads
  select(-denoisedR) %>% 
  pivot_longer(cols = !ID, names_to = "step", values_to = "num_seqs") 

# merge
all <- trkv %>% rbind(both) %>% 
  # order
  mutate(step = factor(step, ordered = TRUE, 
                        levels = c("raw", "adapt", "trim",
                                   "denoisedF", "merged", "nonchim")))
## summarize for supp info 

s <- all %>% group_by(step) %>% get_summary_stats(type = "common")

# write
#write.table(s, file = "data/read-loss-summary.txt", sep = "\t", row.names = FALSE)

# get negative controls
negs <- all %>% filter(ID %in% c("M11", "M148", "M94")) %>% 
  group_by(step) %>% get_summary_stats(type = "common")

# get non-negatives
non <- all %>% filter(!ID %in% c("M11", "M148", "M94")) %>% 
  group_by(step) %>% get_summary_stats(type = "common")

# get positive controls
pos <- all %>% filter(ID %in% c("M10", "M147", "M92", "M93")) %>% 
  group_by(step) %>% get_summary_stats(type = "common")


## ---- visualize ----

ggdensity(all, x = "num_seqs", facet.by = "step")


# get quartiles
all %>% group_by(step) %>% get_summary_stats(num_seqs, type = "quantile")

# get summary stats
all %>% group_by(step) %>% get_summary_stats(num_seqs, type = "mean_sd")

# get sum of raw and nonchim reads
all %>% group_by(step) %>% summarize(sum = sum(num_seqs))

## ---- total lost from raw reads ----

# wrangle
tot <- all %>% 
  filter(step %in% c("raw", "nonchim")) %>% 
  pivot_wider(names_from = "step", values_from = "num_seqs") %>% 
  mutate(total.lost = raw - nonchim,
         perc.lost = ((raw-nonchim)/raw) * 100)


## visualize
ggdensity(tot, x = "perc.lost")
ggdensity(tot, x = "total.lost")

# get summary stats
tot %>% get_summary_stats(perc.lost, type = "quantile")

## ---- where are the most reads lost ----

# make horizontal
calc <- all %>% 
  pivot_wider(names_from = c(step, direction), values_from = num_seqs) %>% 
  ## calculate percent losses
  mutate(ad.loss = ((raw_forward - adapt_forward) / raw_forward) * 100,
         trim.loss = ((adapt_forward - trim_forward) / adapt_forward) * 100,
         denoiseF.loss = ((trim_forward - denoisedF_forward) / trim_forward) * 100,
         denoiseR.loss = ((trim_reverse - denoisedR_reverse) / trim_reverse) * 100,
         merge.loss = ((denoisedF_forward - merged_NA) / denoisedF_forward) * 100,
         chim.loss = ((merged_NA - nonchim_NA) / merged_NA) * 100) %>% 
  # keep only calculated columns
  select(ID, ends_with("loss"))

# make vertical to get summary stats
calc %>% 
  pivot_longer(!ID, names_to = "param", values_to = "percent.loss") %>% 
  group_by(param) %>% 
  get_summary_stats(percent.loss, type = "quantile")
