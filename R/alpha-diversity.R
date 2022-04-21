## alpha diversity analysis

# EVS 12/2021
# 2/2022 model all and account for multiple comparisons
# 3/2022 update: only Simpson (evenness), exploratory for Observed (richness)


require(tidyverse)
require(phyloseq)
require(microViz)
require(rstatix)
require(ggpubr)

load("./data/ps-decontam-filtered-counts.RData")
ps <- pscount 

## ---- get alpha diversity metrics -----

sampdf <- sample_data(ps) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID")

# get diversity
alpha <- estimate_richness(ps, measures = c("Simpson", "Observed")) %>% 
  # make ID
  rownames_to_column(var = "ID") %>% 
  merge(sampdf, by = "ID") 

# get summary statistics 
sum <- alpha %>% 
  group_by(Treatment, Age) %>% get_summary_stats(Simpson, type = "mean_sd")
#write.table(sum, file = "./data/simpson-sum-stats.txt", sep = "\t", row.names = FALSE)

## ----- Observed ----

# exploratory histogram
ggdensity(alpha, x = "Observed", fill = "Treatment",
          facet.by = "Age") # potential outlier in 50 weeks 50 mg/kg
# exploratory boxplot
ggboxplot(alpha, x = "Treatment", y = "Observed", facet.by = "Age")

## ---- Simpson ----

# exploratory plots
hist(alpha$Simpson)
ggboxplot(alpha, x = "Treatment", y = "Simpson", facet.by = "Age")
ggdensity(alpha, x = "Simpson", fill = "Treatment", facet.by = "Age")
# one low potential outlier
which(alpha$Simpson < 0.8)
alpha[64,] # 50 mg/kg at 40 weeks (Bird 221)

# get summary stats to decide on outlier
st <- alpha %>% get_summary_stats(Simpson, type = "mean_sd") %>% 
  mutate(sdx3 = sd * 3,
         outhi = mean + sdx3,
         outlo = mean - sdx3)
# SD of Simpson is 0.03
alpha$Simpson[81] < st$outlo # not low enough to remove


# build test model for linear assumption
testm <- lm(Simpson ~ Treatment, data = alpha)
hist(resid(testm))
qqnorm(resid(testm))
qqline(resid(testm))

## try log transformation
test1 <- lm(log1p(Simpson) ~ Treatment, data = alpha)
hist(resid(test1))
qqnorm(resid(test1))
qqline(resid(test1))

# not a lot of improvement - will stick to un-transformed

##### MODELS #####
## correct for multiple comparisons: 3 hypothesis tests
# get new alpha
alpha.crit <- 0.05 / 3


## test 1: is there a difference between control and T4 (control and metformin)?
mod1 <- t.test(Simpson ~ Treatment, data = filter(alpha, Treatment %in% c("0 mg/kg", "75 mg/kg")))
mod1
# no difference
# adjust p value
p.adjust(mod1$p.value, method = "bonferroni", n = 3)

## test 2: is there a difference in metformin treatment over Age
mod2 <- aov(Simpson ~ Age, data = filter(alpha, !Treatment %in% c("0 mg/kg")))
summary(mod2)
# no difference
# adjust p value
p.adjust(0.121, method = "bonferroni", n = 3)


## sanity check: is there a difference in control over Age
mod3 <- aov(Simpson ~ Age, data = filter(alpha, Treatment == "0 mg/kg"))
summary(mod3) # no difference

## test 3: is there a dose response
mod4 <- aov(Simpson ~ Treatment, data = filter(alpha, Treatment != "0 mg/kg"))
summary(mod4) # none
# adjust p value
p.adjust(0.309, method = "bonferroni", n = 3)
