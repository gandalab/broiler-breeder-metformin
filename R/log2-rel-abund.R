## rel abundance - log2 linear model 

require(microViz)
require(tidyverse)
require(phyloseq)

load("data/ps-decontam-filtered-counts.RData")
ps <- pscount 

## BONFERRONI CORRECTION
## three tests and 805 total genera in all tests
alpha <- 0.05 / 805

## ---- test1: 0 mg/kg v t4 ----

# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("0 mg/kg", "75 mg/kg")) 

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

## adjust ALL tests by hand
## significant after Bonferonni correction?
any(lm_stats$taxatree_stats$p.value < alpha)
# get modified p value
lm_stats$taxatree_stats <- lm_stats$taxatree_stats %>% 
  mutate(p.adj.Bon = 
           p.adjust(p.value, method = "bonferroni", n = 805))

# how many are still significant?
#View(lm_stats$taxatree_stats %>% filter(p.adj.Bon < 0.05))

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)


# save for comparison
mod1 <- lm_stats

## ---- test2: difference between metformin over time ----

# only pairwise comparisons, so use times 1 and 3
# get data
phylo <- ps %>% 
  ps_filter(!Treatment %in% c("0 mg/kg")) %>% 
  ps_filter(Age %in% c("40 weeks", "60 weeks")) %>%
  #ps_mutate(Age = factor(Age, ordered = TRUE, levels = c("40 weeks", "60 weeks"))) %>% 
  tax_fix()

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Age")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

## adjust ALL tests by hand
## significant after Bonferonni correction?
any(lm_stats$taxatree_stats$p.value < alpha)
# get modified p value
lm_stats$taxatree_stats <- lm_stats$taxatree_stats %>% 
  mutate(p.adj.Bon = 
           p.adjust(p.value, method = "bonferroni", n = 805))

# are any still significant?
#View(lm_stats$taxatree_stats %>% filter(p.adj.Bon < 0.05)) # one

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.Bon < 0.05])

# save for later
mod2 <- lm_stats

### ---- sanity check: 0 mg/kg over time ----

# get data
phylo <- ps %>% 
  ps_filter(Treatment %in% "0 mg/kg") %>% 
  ps_filter(Age %in% c("40 weeks", "60 weeks")) %>% 
  tax_fix()

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Age")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

## adjust ALL tests by hand
## significant after Bonferonni correction?
any(lm_stats$taxatree_stats$p.value < alpha) # none
# get modified p value
lm_stats$taxatree_stats <- lm_stats$taxatree_stats %>% 
  mutate(p.adj.Bon = 
           p.adjust(p.value, method = "bonferroni", n = 805))

# save for later
mod2.5 <- lm_stats

## ---- test 3: dose response ----

# get data
phylo <- ps %>% 
  ps_filter(!Treatment %in% c("0 mg/kg")) %>% 
  # pairwise comparisons; 25 mg/kg v 575mg/kg
  ps_filter(!Treatment == "50 mg/kg") %>% 
  tax_fix()

# model on Genus
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform(
    trans = "log2", 
    # chain allows for multiple transformations (TSS then log2)
    chain = TRUE, 
    # for log, replace 0's with half of min value
    zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = "Genus",
    variables = c("Treatment")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

## adjust ALL tests by hand
## significant after Bonferonni correction?
any(lm_stats$taxatree_stats$p.value < alpha)
# get modified p value
lm_stats$taxatree_stats <- lm_stats$taxatree_stats %>% 
  mutate(p.adj.Bon = 
           p.adjust(p.value, method = "bonferroni", n = 805))

# are any still significant?
#View(lm_stats$taxatree_stats %>% filter(p.adj.Bon < 0.05))  # THREE

# what are the effect sizes?
#hist(lm_stats$taxatree_stats$estimate)

# save for later
mod3 <- lm_stats

## ---- total number genera in all comps ----

mod1 # 266
mod2 # 269
mod3 #270

266+269+270 # 805

## ---- save tables for summary stats ----

# model 1
m1 <- mod1$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "75 mg/kg vs 0 mg/kg") %>% ## where POSITIVE is higher in 75 mg/kg
  ungroup() %>% 
  select(-c(taxon, rank, term)) %>% 
  relocate(model, Genus)


# model 2
m2 <- mod2$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "Metformin-treated hens: 60 vs 40 weeks") %>% ## where POSITIVE is higher in 60 weeks
  ungroup() %>% 
  select(-c(taxon, rank, term)) %>% 
  relocate(model, Genus)

# model 2.5
m2.5 <- mod2.5$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "0 mg/kg hens: 60 vs 40 weeks") %>% ## where POSITIVE is higher in 60 weeks
  ungroup() %>% 
  select(-c(taxon,  rank, term)) %>% 
  relocate(model, Genus)


# model 3
m3 <- mod3$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "75 mg/kg vs 25 mg/kg") %>% ## where POSITIVE is higher in 75 mg/kg than 25 mg/kg
  ungroup() %>% 
  select(-c(taxon,  rank, term)) %>% 
  relocate(model, Genus)

# combie and save
all <- rbind(m1, m2, m2.5, m3)

write.table(all, file = "data/rel-abund-results.txt", sep = "\t", row.names = FALSE)
