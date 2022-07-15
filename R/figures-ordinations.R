## Figure 2: PCA on CLR transformation

# EVS 2/2022

require(microViz)
require(colorspace)
require(tidyverse)
require(ggpubr)
require(RColorBrewer)

load("data/ps-decontam-filtered-counts.RData")
ps <- pscount %>% 
  ps_mutate(Metformin = if_else(Treatment == "0 mg/kg", "Control", "Metformin"))

## get colors
#treatcols <- qualitative_hcl("Set3", n = 4)
treatcols <- brewer.pal(n = 4, name = "Dark2")
names(treatcols) <- levels(samdat_tbl(ps)$Treatment)
#agecols <- qualitative_hcl("Dynamic", n = 3)
agecols <- brewer.pal(n = 3, name = "Set1")
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
           tax_lab_style = tax_lab_style(size = 4,
                                         justify = "side",
                                         type = "label",
                                         position = position_jitter(height = 0.7)),
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
ggsave(filename = "R/figures/pca-control-75mg-ordination.png", plot = last_plot(), dpi = 600)

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
           tax_lab_style = tax_lab_style(size = 3,
                                         justify = "side"),
           auto_caption = NA) + 
  stat_ellipse(aes(color = Age), alpha = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_color_manual(values = agecols) +
  theme(
    text = element_text(size = 16))

# save
ggsave(filename = "R/figures/pca-age.png", plot = last_plot(), dpi = 600)

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
           tax_lab_style = tax_lab_style(size = 4,
                                         justify = "side",
                                         position = position_jitter(height = 0.8,
                                                                    width = 0.2)),
           auto_caption = NA) + 
  stat_ellipse(aes(color = Treatment), alpha = 0.5) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_color_manual(values = c(treatcols[2], treatcols[3], treatcols[4])) +
  theme(
    text = element_text(size = 16))

# save
ggsave(filename = "R/figures/pca-dose.png", plot = last_plot(), dpi = 600)


## ---- supp figure 1: beta dispersions ----

## dose response

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
png(filename = "R/figures/SuppFig_doseresp_betadisp.png")
# plot
plot(bd$Treatment$model, main = NULL, sub = NULL)
title("All metformin doses (25, 50, and 75 mg/kg", adj = 0)
dev.off()
