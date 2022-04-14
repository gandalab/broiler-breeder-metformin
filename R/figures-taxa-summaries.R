## FIGURE 1: taxa summary

# EVS 2/2022

require(microViz)
require(phyloseq)
require(ggpubr)
require(tidyverse)
require(RColorBrewer)
require(colorspace)

load("data/ps-decontam-filtered-counts.RData")


ps <- pscount

# assign colors
cols <- sequential_hcl(palette = "TealGrn", n = 4)
#names(cols) <- unique(samdat_tbl(ps)$Treatment)
names(cols) <- c("75 mg/kg", "50 mg/kg", "25 mg/kg", "0 mg/kg")

# assign diverging colors
cols2 <- qualitative_hcl(palette = "Dark3", n = 2)
names(cols2) <- c("Control", "Metformin")

# get RColorBrewer - Set2
colsSet2 <- brewer.pal(n = 3, name = "Set2")[1:2]
names(colsSet2) <- c("Control", "Metformin")

## ---- barplot ----
## this has been moved to Supp data (Supp Figure 2)

ps %>% 
  tax_fix() %>% 
  comp_barplot(
    tax_level = "Phylum",
    label = NULL, # don't label each sample (too many!)
    n_taxa = length(get_taxa_unique(ps, "Phylum")), # plot all taxa
    bar_outline_color = NA,
    bar_width = 0.98
  ) +
  theme(text = element_text(size = 16)) +
  labs(y = "Relative Abundance")

# save
ggsave(filename = "R/figures/allsamps-barplot-phylum.png", plot = last_plot(), dpi = 600)

## ---- top genera ----

# plot top 10 most abundant genera

top <- prune_taxa(names(sort(taxa_sums(ps), TRUE))[1:10], ps)

# (this gives a bunch of warnings due to merging)
p <- top %>% 
  merge_samples(group = "Treatment") %>% 
  comp_barplot(tax_level = "Genus",
               n_taxa = 10,
               bar_width = 0.95,
               bar_outline_colour = "grey5",
               sample_order = c("0 mg/kg", "25 mg/kg", "50 mg/kg", "75 mg/kg")
  ) +
  theme_pubr(legend = "right", x.text.angle = 45) +
  labs(y = "Relative Abundance") +
  theme(text = element_text(size = 16))

## save with normal legend
ggsave(filename = "R/figures/genus-treatment-relabun.png", plot = p, dpi = 600)

## which ranks do these belong to?
tax_table(prune_taxa(names(sort(taxa_sums(psg), TRUE))[1:6], psg))

## ---- heatmap ----

# recode metformin for coloring 
ps <- ps %>% 
  # recode metformin treatment
  ps_mutate(
            # combine metformin treatments
            Metformin = if_else(Treatment == "0 mg/kg", "Control", "Metformin")) 


## for heatmap: view color palettes from colorspace: https://colorspace.r-forge.r-project.org/articles/colorspace.html

p <- ps %>% 
  tax_transform("clr", rank = "Phylum") %>% 
  #tax_filter(min_prevalence = 0.5) %>%  # keep only taxa in most samples 
  comp_heatmap(
    colors = heat_palette(palette = "Green-Brown", rev = TRUE), name = "CLR",
    tax_anno = taxAnnotation(
      Prevalence = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
    )
    # take out treatment label,
    #sample_anno = sampleAnnotation(
    # Treatment = anno_sample("Metformin"),
    #col = list(Metformin = colsSet2), border = FALSE
    #                      )
  )

# save - complex heatmap so save as pdf
pdf(file = "R/figures/heatmap-phyla.pdf")
p
dev.off()
