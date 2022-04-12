## rel abundance - log2 linear model 

require(microViz)
require(tidyverse)
require(phyloseq)

load("data/ps-decontam-filtered-counts.RData")
ps <- pscount %>% 
  ps_mutate(Time = as.factor(recode(Sample.Date, '8/8/19' = '40 weeks', '10/16/19' = '50 weeks', '12/20/19' = '60 weeks'))) %>% 
  tax_fix()

## BONFERRONI CORRECTION
## three tests and 789 total genera in all tests
alpha <- 0.05 / (789)

## ---- test1: control v t4 ----

# subset 
phylo <- ps %>% 
  ps_filter(Treatment %in% c("Control", "T4")) 

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

# adjust for multiple comparisons (within the test)
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

## adjust ALL tests by hand
## significant after Bonferonni correction?
any(lm_stats$taxatree_stats$p.value < alpha)
# get modified p value
lm_stats$taxatree_stats <- lm_stats$taxatree_stats %>% 
  mutate(p.adj.Bon = 
              p.adjust(p.value, method = "bonferroni", n = 790))

# how many are still significant?
View(lm_stats$taxatree_stats %>% filter(p.adj.Bon < 0.05))

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.Bon < 0.05])

# save for comparison
mod1 <- lm_stats

## ---- test2: difference between metformin over time ----

# only pairwise comparisons, so use times 1 and 3
# get data
phylo <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
  ps_filter(Time %in% c("40 weeks", "60 weeks")) %>%
  ps_mutate(Age = factor(Time, ordered = TRUE, levels = c("40 weeks", "60 weeks"))) %>% 
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
    variables = c("Time")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

## adjust ALL tests by hand
## significant after Bonferonni correction?
any(lm_stats$taxatree_stats$p.value < alpha)
# get modified p value
lm_stats$taxatree_stats <- lm_stats$taxatree_stats %>% 
  mutate(p.adj.Bon = 
           p.adjust(p.value, method = "bonferroni", n = 790))

# are any still significant?
View(lm_stats$taxatree_stats %>% filter(p.adj.Bon < 0.05)) # none

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for later
mod2 <- lm_stats

### ---- sanity check: control over time ----

# get data
phylo <- ps %>% 
  ps_filter(Treatment %in% "Control") %>% 
  ps_filter(Time %in% c("40 weeks", "60 weeks")) %>% 
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
    variables = c("Time")
  )

# get stats
lm_stats <- taxatree_models2stats(lm_models)
lm_stats$taxatree_stats

# adjust for multiple comparisons
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

## adjust ALL tests by hand
## significant after Bonferonni correction?
any(lm_stats$taxatree_stats$p.value < alpha) # none
# get modified p value
lm_stats$taxatree_stats <- lm_stats$taxatree_stats %>% 
  mutate(p.adj.Bon = 
           p.adjust(p.value, method = "bonferroni", n = 801))

# are any still significant?
lm_stats$taxatree_stats %>% filter(p.adj.Bon < 0.05) # none

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for later
mod2.5 <- lm_stats

## ---- test 3: dose response ----

# get data
phylo <- ps %>% 
  ps_filter(!Treatment %in% c("Control", "Males")) %>% 
  # pairwise comparisons; T2 v T4
  ps_filter(!Treatment == "T3") %>% 
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

# adjust for multiple comparisons
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

# which are significant?
lm_stats$taxatree_stats %>% filter(p.adj.BH.rank < 0.05)

## adjust ALL tests by hand
## significant after Bonferonni correction?
any(lm_stats$taxatree_stats$p.value < alpha)
# get modified p value
lm_stats$taxatree_stats <- lm_stats$taxatree_stats %>% 
  mutate(p.adj.Bon = 
           p.adjust(p.value, method = "bonferroni", n = 790))

# are any still significant?
View(lm_stats$taxatree_stats %>% filter(p.adj.Bon < 0.05))  # ONE

# what are the effect sizes?
hist(lm_stats$taxatree_stats$estimate)

# which taxa are changing?
unique(lm_stats$taxatree_stats$taxon[lm_stats$taxatree_stats$p.adj.BH.rank < 0.05])

# save for later
mod3 <- lm_stats

## ---- total number genera in all comps ----

mod1 # 262
mod2 # 264
mod3 #263

## ---- save tables for summary stats ----

# model 1
m1 <- mod1$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "75 mg/kg versus 0 mg/kg") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m1, file = "data/rel-abund-control_T4-results.txt", sep = "\t", row.names = FALSE)

# model 2
m2 <- mod2$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "Metformin-treated hens: Age") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m2, file = "data/rel-abund-age-results.txt", sep = "\t", row.names = FALSE)

# model 2.5
m2.5 <- mod2.5$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "Control hens: Age") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m2.5, file = "data/rel-abund-Age-control-results.txt", sep = "\t", row.names = FALSE)

# model 3
m3 <- mod3$taxatree_stats %>% 
  mutate(Genus = str_remove(taxon, "G: "),
         model = "Dose response") %>% 
  ungroup() %>% 
  select(-c(taxon, p.adj.BH.rank, rank, term)) %>% 
  relocate(model, Genus)

write.table(m3, file = "data/rel-abund-dose-results.txt", sep = "\t", row.names = FALSE)

