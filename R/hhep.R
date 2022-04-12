## HHEP data
#EVS 3/2022

library(tidyverse)
library(rstatix)
library(microViz)
require(ggpubr)
require(colorspace)

# get colors
cols <- qualitative_hcl("Dark2", n = 4)


# get egg data
df <- readxl::read_xlsx(path = "../egg-production-for-R.xlsx",
                        sheet = "hhep_compiled") %>% 
  # make vertical
  pivot_longer(cols = !age_weeks, names_to = "Treatment", values_to = "hhep")

## diagnostics: 
# fit linear model
testmod <- lm(hhep ~ age_weeks * Treatment, data = df)

## lump into times
df <- df %>% 
  mutate(period = case_when(
    age_weeks < 25 ~ "0",
    age_weeks >= 25 & age_weeks < 40 ~ "1",
    age_weeks >= 40 & age_weeks < 50 ~ "2",
    age_weeks >= 50 & age_weeks < 70 ~ "3"
  ))

## plot
ggdensity(df, x = "hhep", color = "period")

## get summary stats
df %>% group_by(Treatment, period) %>% get_summary_stats(hhep, type = "common")

# remove time 0 - few data points 
dat <- df %>% filter(!period == 0)

# test model
testmod <- lm(hhep ~ Treatment * period, data = dat)

# run model diagnostics
hist(resid(testmod)) 
qqnorm(resid(testmod))
qqline(resid(testmod))

### ---- anova -----

## mirroring Evelyn's analysis, one test for each time point

## period 1
summary(aov(hhep ~ Treatment, data = filter(dat, period == 1))) # no effect


## period 2
summary(aov(hhep ~ Treatment, data = filter(dat, period == 2))) # no effect

## period 3
summary(aov(hhep ~ Treatment, data = filter(dat, period == 3)))# effect
mod <- aov(hhep ~ Treatment, data = filter(dat, period == 3))

# Tukey's post hoc
TukeyHSD(mod) # Control -T3 only 

### ---- plot ----

# rename treatments
plot.dat <- df %>% 
  mutate(treat = recode_factor(Treatment,
                               Control = "0 mg/kg",
                               T2 = "25 mg/kg",
                               T3 = "50 mg/kg",
                               T4 = "75 mg/kg",
                               .ordered = TRUE),
         age_weeks = as.integer(age_weeks))

p <- ggline(plot.dat, y = "hhep", x = "age_weeks", color = "treat",
       plot_type = "l", # line only
       size = 1,
       xlab = "Age (weeks)", ylab = "HHEP",
       numeric.x.axis = TRUE
       ) +
  
  scale_y_continuous(expand = expansion(add = 15)) +
  labs(color = "") +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = 25, linetype = 3) +
  geom_vline(xintercept = 40, linetype = 2) +
  geom_vline(xintercept = 50, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  geom_label(x = 29, y = 102, label = "metformin starts") +
  geom_label(x = 42, y = 102, label = "1st sample") +
  geom_label(x = 52, y = 102, label = "2nd sample") +
  geom_label(x = 62, y = 102, label = "3rd sample")

ggpar(p, x.text.angle = 45, xticks.by = 5)

# save
ggsave(filename = "R/figures/hhep.png", plot = last_plot(), dpi = 600)


## ---- summary stats ----

# how many hens were laying in each treatment at 65 weeks
df %>% filter(age_weeks == 65) 

# 40, 50, and 60 weeks (sampling times)
df %>% filter(age_weeks %in% c(40, 50, 60))
