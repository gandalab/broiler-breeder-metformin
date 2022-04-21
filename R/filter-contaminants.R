### contaminants and filtering
# EVS UPDATE 4/22 - nuke based on pos controls, not decontam

library(tidyverse)
library(phyloseq)
library(microViz)
library(ggpubr)
library(vegan)

## get raw phyloseq
load("data/raw-ps.RData")

## REMOVE: Males and DNA pos control (miscommunication with Novogene)
ps <- ps %>% 
  ps_filter(!Treatment == "Males") %>% 
  ps_filter(!Treatment == "PCISO1")

## how many taxa in the negative controls?
neg <- ps %>% ps_filter(str_detect(Treatment, "NC")) # 1577
pos <- ps %>% ps_filter(str_detect(Treatment, "PC"))

## ---- build positive Zymo control ----

## control used in DNA extraction: ZymoBiomics Community Standard Lot ZRC190633, Cat D6300

# strain info for 16S here: https://files.zymoresearch.com/pdf/d6300-_zymobiomics_microbial_community_standard_v1-1-3.pdf

# build taxonomy -- some names have a slight change to match Silva taxonomy
tax <- data.frame(
  #row.names = c("Kingdom", "Genus"),
  p.aeru = c("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas"),
  e.coli = c("Bacteria","Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia-Shigella"),
  s.ent = c("Bacteria","Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Salmonella"),
  l.ferm = c("Bacteria","Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"),
  e.fae = c("Bacteria","Firmicutes", "Bacilli", "Lactobacillales", "Enterococcaceae", "Enterococcus"),
  s.aur = c("Bacteria","Firmicutes", "Bacilli", "Staphylococcales", "Staphylococcaceae", "Staphylococcus"),## CHANGED THE ORDER
  l.mon = c("Bacteria","Firmicutes", "Bacilli", "Lactobacillales", "Listeriaceae", "Listeria"), ## CHANGED THE ORDER
  b.sub = c("Bacteria","Firmicutes", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus")
  
)  %>% t()

colnames(tax) <-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# build Zymo mock community by hand
zymo <- phyloseq(otu_table(
  data.frame(
    row.names = "ZymoComm",
    p.aeru = 0.042,
    e.coli = 0.101,
    s.ent = 0.104,
    l.ferm = 0.184,
    e.fae = 0.099,
    s.aur = 0.155,
    l.mon = 0.141,
    b.sub = 0.174
  ), taxa_are_rows = FALSE
),
tax_table(
  tax
  
),
# have to do this for microViz...
sample_data(data.frame(
  row.names = "ZymoComm",
  Treatment = "Zymo"
))
)

## ---- plot pos and neg controls with Zymo ----

# get zymo
zym <- ps_melt(zymo) %>% 
  mutate(all.name = paste(Phylum, Class, Order, Family, Genus, sep = "_")) 

## get raw data and add Zymo
raw <- ps %>% 
  tax_glom("Genus") %>% 
  ps_melt() %>% 
  mutate(all.name = paste(Phylum, Class, Order, Family, Genus, sep = "_")) %>% 
  full_join(zym)

## plot
ggplot(data = filter(raw, Treatment %in% c("NCP1", "NCP2","PCP1", "PCP2", "Zymo")), mapping = aes_string(x = "Treatment", y = "Abundance")) +
  geom_bar(aes(fill = Class), stat = "identity", position = "fill") +
  labs(x = NULL, y = "Relative Abundance")

## ---- filter contaminants ----

## REMOVE: ALL TAXA IN NEG CONTROLS AND NON-MOCK COMM TAXA IN POS CONTROLS

# get list of mock community taxa
mockcomm <- data.frame(tax) %>% 
  mutate(all.name = paste(Phylum, Class, Order, Family, Genus, sep = "_")) %>% 
  select(all.name) 
mock <- mockcomm$all.name

## first, detect whether the taxa belongs in the mock community
ismock <- raw %>% 
  mutate(in.mock = if_else(all.name %in% mock, TRUE, FALSE)) %>% 
  mutate(is.control = if_else(Treatment %in% c("Males", "Control", "T2", "T3", "T4"), FALSE, TRUE)) %>% 
  mutate(is.pos = if_else(str_detect(Treatment, "PCP"), TRUE, FALSE),
         is.neg = if_else(str_detect(Treatment, "NCP"), TRUE, FALSE),
         is.zym = if_else(Treatment == "Zymo", TRUE, FALSE))

# filter based on all controls and mock community
filt <- ismock %>% 
  mutate(is.contam = case_when(
    is.neg ~ TRUE,
    !is.neg & !is.control ~ FALSE,
    !is.neg & is.zym ~ FALSE,
    !is.neg & in.mock & is.pos ~ FALSE,
    !is.neg & !in.mock & is.pos ~ TRUE
  )) 

# 3,114 "contaminants"
length(which(filt$is.contam))
length(which(ismock$is.neg)) # 1,569 in neg control

### replicate decontam plot - what is prevalence of 'contaminants' in controls vs samples ###

# convert to presence/absence to get prevalence 
tsums <- filt %>% 
  mutate(pres = if_else(Abundance > 0, 1, 0)) %>% 
  group_by(is.control, OTU) %>% 
  summarize(taxa.prev = sum(pres))
# get contaminants 
contams <- filt %>% select(OTU, is.contam) %>% distinct()
df.pa <- tsums %>% 
  left_join(contams) %>% 
  pivot_wider(names_from = "is.control", values_from = "taxa.prev") %>% 
  rename(prev.controls = `TRUE`, 
         prev.samps = `FALSE`) %>% 
  drop_na() # artifact of Zymo strains
ggplot(data = df.pa, aes(x = prev.controls, y = prev.samps, color = is.contam)) +
  geom_point() +
  labs(x = "Prevalence in controls", y = "Prevalence in samples")

# remove contaminants
decon <- filt %>% 
  filter(!is.contam)

## ---- filter rel abun, etc ----

# filter based on relative abundance, prevalence, and bacteria
decon1 <- decon %>% 
  filter(Kingdom == "Bacteria") %>% 
  filter(!is.na(Phylum))

#### get data back to phyloseq and filter for relative abundance
newotu <- otu_table(decon1 %>% 
                      select(Sample, OTU, Abundance) %>% 
                      pivot_wider(names_from = OTU, values_from = Abundance, values_fill = 0) %>% 
                      column_to_rownames(var = "Sample"), taxa_are_rows = FALSE)
newtax <- decon1 %>% 
  select(OTU, Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  distinct() %>% 
  column_to_rownames(var = "OTU") %>% 
  as.matrix() 
newtax <- tax_table(newtax)
newsamp <- decon1 %>% 
  select(Sample, Bird, Treatment, Sample.Date, is.control) %>% 
  distinct() %>% 
  column_to_rownames(var = "Sample")

newps <- phyloseq(newotu, newtax,
                  sample_data(newsamp))

## filter for relative abundance
newps <- tax_fix(newps)

# transform to relative abundance
pst <- transform_sample_counts(newps, function(x) x / sum(x) )

# remove taxa with total relative abundance less than 10e-5
psr <- filter_taxa(pst, function(x) mean(x) > 1e-5, TRUE)

## get count data from raw phyloseq and subset to contaminants
# this will remove the 8 strains from Zymo control 
pscount <- subset_taxa(ps, taxa_names(ps) %in% taxa_names(psr))
# remove sequencing controls and fix treatment labels for downstream analysis
pscount <- pscount %>% 
  ps_filter(!Bird == "Sequence-Control") %>% 
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
                                      .ordered = TRUE))

# save phyloseq
#save(pscount, file = "data/ps-decontam-filtered-counts.RData")
