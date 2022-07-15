### Figures for poster 

# EVS 5/2022

### ---- relative abundnce in 3 panels ----

require(microViz)
require(tidyverse)
require(phyloseq)
require(ggpubr)
require(colorspace)

load("data/ps-decontam-filtered-counts.RData")
ps <- pscount

## get model info
mods <- read.table("data/rel-abund-results.txt", sep = "\t", header = TRUE)

# get colors
cols <- qualitative_hcl("Dark2", n = 3)


# get significant taxa for each model
df <- mods %>% 
  filter(p.adj.Bon < 0.05) %>% 
  # shorten name of model
  mutate(model = case_when(
    str_detect(model, "Metformin-treated") ~ "40 vs 60 weeks",
    str_detect(model, "0 mg/kg") ~ "0 vs 75 mg/kg",
    str_detect(model, "25 mg/kg") ~ "25 vs 75 mg/kg"
  ))

# get genera
taxa <- unique(df$Genus)

## get full taxonomic names for the ignificant genera
names <- ps %>% tax_select(taxa, ranks_searched = "Genus", n_typos = 0) %>% 
  ps_melt() %>% 
  select(Phylum, Class, Order, Family, Genus) %>% 
  distinct() %>% 
  # make name from Class & Genus
  mutate(taxa = paste(Family, Genus, sep = " "))

## join
plot_dat <- df %>% full_join(names) %>% 
  ## FLIP THE SIGNS TO BE MORE INTUITIVE 
  mutate(t.estimate = -estimate) %>% 
  # fix double name
  mutate(taxa = if_else(str_detect(taxa, "UCG"), "Lachnospiraceae UCG-010", taxa)) 

### FACET WRAP
p <- ggdotchart(plot_dat, x = "taxa", y = "t.estimate", group = "model",
           facet.by = "model", 
                color = "model", shape = "model",
                sorting = "descending",
                add = "segments", add.params = list(size = 4, alpha = 0.5),
                dot.size = 7,
                # add log2 fold change in the dot
                label = round(plot_dat$t.estimate, 2),
                font.label = list(color = "black", size = 18, vjust = -1.5),
                rotate = TRUE,
                xlab = "",
                ylab = "log2 fold change",
           legend = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # remember this is flipped, so scale_y is actually the x axis
  scale_y_continuous(expand = expansion(mult = c(.1, .15))) +
  # remove legend label
  #labs(color = "", shape = "")  +
  scale_color_manual(values = cols) +
  theme(text = element_text(size = 20))

# change font
ggpar(p, font.ytickslab = "italic", orientation = "vertical") +
  #guides(color = guide_legend(title.hjust = -3)) +
  scale_x_discrete(position = "top") +
  theme(strip.background =element_rect(fill="white"))

# save
ggsave(filename = "R/figures/poster-rel-abund.png", plot = last_plot(), dpi = 600,
       height = 6, width = 11, units = "in")

