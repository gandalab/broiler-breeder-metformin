## beta diversity

# Emily Van Syoc 12/2021
# 2/2022 - replace bray with aitchison and Bonferroni correction

require(tidyverse)
require(phyloseq)
require(ggpubr)
require(vegan)
require(microViz)
require(pairwiseAdonis)

load("./data/ps-decontam-filtered-counts.RData")
ps <- pscount

### THREE MODELS: BONFERRONI CORRECTION
# get new critical alpha value
alpha <- 0.05 / 3 

## --- global test for interactions ----

# EXPLORATORY ONLY
ps %>% 
  tax_fix() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  dist_calc("aitchison") %>% 
  dist_permanova(
    variables = c("Treatment", "Age"),
    interactions = "Treatment*Age",
    perm = 9999
  )

### ---- test 1: difference between 0 mg/kg and T4 ----

# get data
t1 <- ps %>% 
  ps_filter(Treatment %in% c("0 mg/kg", "75 mg/kg")) %>% 
  # fix ordered factor
  ps_mutate(Treatment = as.character(Treatment)) %>% 

  # transform to relative abundance
  tax_transform("compositional", rank = "Genus", keep_counts = FALSE) %>% 
  #dist_calc("bray") %>% 
  dist_calc("aitchison")

# test beta dispersion
bd <- t1 %>% dist_bdisp(variables = "Treatment") %>% bdisp_get() # no significance
# plot
plot(bd$Treatment$model, main = NULL)

# test
mod1 <- t1 %>% 
  dist_permanova(
    seed = 123,
    variables = c("Treatment"),
    n_perms = 9999
  ) # significant

## significant after Bonferonni correction?
mod1$permanova$`Pr(>F)`[1] < alpha # TRUE
# get modified p value
p.adjust(mod1$permanova$`Pr(>F)`[1], method = "bonferroni", n = 3)

# plot
mod1 %>% 
  ord_calc(method = "auto") %>% 
  ord_plot(color = "Treatment") +
  stat_ellipse(aes(linetype = Treatment, color = Treatment)) +
  #geom_text(aes(x = -0.8, y = -1.4), label = "adonis p value = 0.0026") +
  ggtitle("75 mg/kg metformin versus 0 mg/kg")

## ---- test 2: difference in metformin treatment over time ----

# get data
t2 <- ps %>% 
  ps_filter(!Treatment %in% c("0 mg/kg")) %>% 
  
  # fix ordered factor
  ps_mutate(Age = as.character(Age)) %>% 
  
  tax_fix() %>% 
  # transform to relative abundance
  tax_transform("compositional", rank = "Genus", keep_counts = FALSE) %>% 
  #dist_calc("bray")
  dist_calc("aitchison")

# test beta dispersion
bd <- t2 %>% dist_bdisp(variables = "Age") %>% bdisp_get() # not significant
# plot
#plot(bd$Sample.Date$model)

# test
mod2 <- t2 %>% 
  dist_permanova(
    seed = 123,
    variables = c("Age"), # not significant
    n_perms = 9999
  )

## significant after Bonferonni correction?
mod2$permanova$`Pr(>F)`[1] < alpha # FALSE
# get modified p value
p.adjust(mod2$permanova$`Pr(>F)`[1], method = "bonferroni", n = 3)

## ---- test 3: metformin treatment dose response ----

# get data
t3 <- ps %>% 
  ps_filter(!Treatment %in% c("0 mg/kg")) %>%
  # fix ordered factor
  ps_mutate(Treatment = as.character(Treatment)) %>% 
  # transform to relative abundance
  tax_transform("compositional", rank = "Genus", keep_counts = FALSE) %>% 
  #dist_calc("bray")
  dist_calc("aitchison")

# test beta dispersion
bd <- t3 %>% dist_bdisp(variables = "Treatment") %>% bdisp_get() # significance
# plot
plot(bd$Treatment$model)

# test
mod3 <- t3 %>% 
  dist_permanova(
    seed = 123,
    variables = c("Treatment"), # significant
    n_perms = 9999
  )

## significant after Bonferonni correction?
mod3$permanova$`Pr(>F)`[1] < alpha # TRUE
# get modified p value
p.adjust(mod3$permanova$`Pr(>F)`[1], method = "bonferroni", n = 3)

# get pairwise comparisons
dis <- dist_get(t3)
sampdf <- samdat_tbl(t3)
pairwise.adonis2(dis ~ Treatment, data = sampdf)


# plot
mod3 %>% 
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "Treatment") +
  stat_ellipse(aes(linetype = Treatment, color = Treatment)) +
  #geom_text(aes(x = -0.8, y = -1.4), label = "adonis p value = 2e-4") +
  ggtitle("Dose response")

