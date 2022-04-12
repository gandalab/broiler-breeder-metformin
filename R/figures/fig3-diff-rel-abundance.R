## Figure 3: Differential relative abundance

# EVS 2/2022

require(microViz)
require(tidyverse)
require(phyloseq)
require(ggpubr)

load("data/ps-decontam-filtered-counts.RData")

# load the model information
source("R/relative-abundance/log2-linearmod.R")

# recode for better names
ps <- pscount %>% 
  ps_mutate(Time = as.factor(recode(Sample.Date, '8/8/19' = '40 weeks', '10/16/19' = '50 weeks', '12/20/19' = '60 weeks')),
            
            Treatment = recode_factor(Treatment,
                                      Control = "0 mg/kg",
                                      T2 = "25 mg/kg",
                                      T3 = "50 mg/kg",
                                      T4 = "75 mg/kg",
                                      .ordered = TRUE)) %>% 
  
  tax_fix()

# get colors
cols <- qualitative_hcl("Dark2", n = 3)

## ---- lollipop plot ----

# combine mod1, mod2, and mod3 (very few significant taxa)
df<- mod1$taxatree_stats %>% 
  # add identifying column
  mutate(model = "75 vs 0 mg/kg") %>% 
  rbind(mod2$taxatree_stats %>% mutate(model = "60 vs 40 weeks")) %>%  
  rbind(mod3$taxatree_stats %>% mutate(model = "75 vs 25 mg/kg")) %>% 
  # get only sig genera
  filter(p.adj.Bon < 0.05)

## get full taxonomic names for the 4 significant genera
names <- mod1 %>% ps_get() %>% tax_select(c("G: Acinetobacter",
                                            "G: Alloiococcus",
                                            "G: Cellulosilyticum",
                                            "G: Lachnospiraceae UCG-010",
                                            "G: [Ruminococcus] gnavus group"),
                                          ranks_searched = "Genus", n_typos = 0) %>% 
  ps_melt() %>% 
  select(Phylum, Class, Order, Family, Genus) %>% 
  distinct() %>% 
  # make name from Class & Genus
  mutate(taxa = paste(
    str_remove(Class, "C: "),
    str_remove(Genus, "G: "),
    sep = " "
  ))

## join
plot_dat <- df %>% full_join(names, by = c("taxon" = "Genus"))

p <- ggdotchart(plot_dat, x = "taxa", y = "estimate",
           color = "model", shape = "model",
           sorting = "descending",
           add = "segments", add.params = list(size = 3, alpha = 0.5),
           dot.size = 5,
           # add log2 fold change in the dot
           label = round(plot_dat$estimate, 2),
           font.label = list(color = "black", size = 11, vjust = -1.5),
           
           rotate = TRUE,
           xlab = "",
           ylab = "log2 fold change") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "", shape = "")  +
  scale_color_manual(values = cols)

# change font
p1 <- ggpar(p, font.ytickslab = "italic") +
  guides(color = guide_legend(title.hjust = -3))

# save
ggsave(filename = "R/figures/fig3-rel-abund.png", plot = p1, dpi = 600)

