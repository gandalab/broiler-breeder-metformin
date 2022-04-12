## Figure 2: PCoA of Bray-Curtis distance 

# EVS 2/2022

require(microViz)
require(colorspace)
require(tidyverse)

load("data/ps-decontam-filtered-counts.RData")
ps <- pscount %>% 
  ps_filter(!Treatment == "Males") %>% 
  ps_mutate(Age = recode_factor(Sample.Date, 
                                '8/8/19' = "40 weeks",
                                '10/16/19' = "50 weeks",
                                '12/20/19' = "60 weeks",
                                .ordered = TRUE),
            Treatment = recode_factor(Treatment,
                                      Control = "0 mg/kg",
                                      T2 = "25 mg/kg",
                                      T3 = "50 mg/kg",
                                      T4 = "75 mg/kg",
                                      .ordered = TRUE),
            # recode binary treatment
            Metformin = if_else(Treatment == "0 mg/kg", "Control", "Metformin")) %>% 
  
  tax_fix()

## get colors
treatcols <- qualitative_hcl("Set3", n = 4)
names(treatcols) <- levels(samdat_tbl(ps)$Treatment)
agecols <- qualitative_hcl("Dynamic", n = 3, l = 70)
names(agecols) <- levels(samdat_tbl(ps)$Age)

## ---- panel A: control v T4 ----

ps %>% 
  ps_filter(Treatment %in% c("0 mg/kg", "75 mg/kg")) %>% 
  # transform to relative abundance
  tax_transform("clr", rank = "Family", keep_counts = FALSE) %>% 
  #dist_calc("bray") %>% 
  # ordinate
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Treatment",
           size = 3,
           # plot taxa loadings 
           plot_taxa = 1:5,
           tax_lab_style = tax_lab_style(size = 3),
           auto_caption = NA) + 
  # add ellipses
  stat_ellipse(aes(color = Treatment), alpha = 0.5) +
  # make plot bigger
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  # color points
  scale_color_manual(values = c(treatcols[1], treatcols[4])) +
  # make text bigger
  theme(
    text = element_text(size = 16))

# save
ggsave(filename = "R/figures/fig2-control-75mg-ordination.png", plot = last_plot(), dpi = 600)
  
## ---- panel B: metformin over time ----


# 2/26 update: add control samples for comparison
ps %>% 
  #ps_filter(!Treatment %in% c("0 mg/kg")) %>% 
  # transform to relative abundance
  tax_transform("clr", rank = "Family", keep_counts = FALSE) %>% 
  #dist_calc("bray") %>% 
  # ordinate
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Age", shape = "Metformin",
           size = 3,
           plot_taxa = 1:5,
           tax_lab_style = tax_lab_style(size = 3),
           auto_caption = NA) + 
  stat_ellipse(aes(color = Age), alpha = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_color_manual(values = agecols) +
  theme(
    text = element_text(size = 16))

# save
ggsave(filename = "R/figures/fig2-age-ordination.png", plot = last_plot(), dpi = 600)

## ---- panel C: dose response ----

ps %>% 
  ps_filter(!Treatment %in% c("0 mg/kg")) %>% 
  # transform to relative abundance
  tax_transform("clr", rank = "Family", keep_counts = FALSE) %>% 
  #dist_calc("bray") %>% 
  # ordinate
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Treatment",
           size = 3,
           plot_taxa = 1:5,
           tax_lab_style = tax_lab_style(size = 3),
           auto_caption = NA) + 
  stat_ellipse(aes(color = Treatment), alpha = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_color_manual(values = c(treatcols[2], treatcols[3], treatcols[4])) +
  theme(
    text = element_text(size = 16))

# save
ggsave(filename = "R/figures/fig2-dose-ordination.png", plot = last_plot(), dpi = 600)


## ---- supp figure 1: beta dispersions ----

## panel A: model 1 (control v t4)

# get data
t1 <- ps %>% 
  ps_mutate(Treatment = as.character(Treatment)) %>% 
  ps_filter(Treatment %in% c("0 mg/kg", "75 mg/kg")) %>% 
  # transform to relative abundance
  tax_transform("compositional", rank = "Genus", keep_counts = FALSE) %>% 
  #dist_calc("bray") %>% 
  dist_calc("aitchison") 
  # fix 
  
  # test beta dispersion
bd <- t1 %>% dist_bdisp(variables = "Treatment") %>% bdisp_get() # significance

# open graphics device
png(filename = "R/figures/SuppFig4_control-t4_betadisp.png")
# plot
plot(bd$Treatment$model, sub = NULL, main = NULL)
title("A. 0 vs 75 mg/kg metformin", adj = 0)
dev.off()

## Panel B: dose response

# get data
t3 <- ps %>% 
  ps_mutate(Treatment = as.character(Treatment)) %>% 
  ps_filter(!Treatment %in% c("0 mg/kg")) %>% 
  # transform to relative abundance
  tax_transform("compositional", rank = "Genus", keep_counts = FALSE) %>% 
  dist_calc("aitchison")

# test beta dispersion
bd <- t3 %>% dist_bdisp(variables = "Treatment") %>% bdisp_get() # significance

# open graphics device
png(filename = "R/figures/SuppFig4_doseresp_betadisp.png")
# plot
plot(bd$Treatment$model, main = NULL, sub = NULL)
title("B. All metformin doses (25, 50, and 75 mg/kg", adj = 0)
dev.off()
