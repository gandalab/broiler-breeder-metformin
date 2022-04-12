# for SRA submission

# EVS 4/2022

library(tidyverse)
library(microViz)

# get sample data from phyloseq object - including controls
load("data/raw-ps.RData")

# remove males (not used in manuscript)
psf <- ps %>% 
  ps_filter(!Treatment == "Males") %>% 
  # remove DNA extraction control
  ps_filter(!Treatment == "PCISO1")

# get sample data
sampdf <- samdat_tbl(psf)

## change Sample.Date to 1/1/2020 for sequencing controls
sampdf <- sampdf %>% 
  mutate(Sample.Date = if_else(Bird == "Sequence-Control", "1/1/2020", Sample.Date)) %>% 
  # change Bird to be unique
  mutate(Bird = if_else(Bird == "Sequence-Control", Treatment, Bird)) %>% 
  # update treatments to be readable
  mutate(Treatment = case_when(
    Treatment %in% "PCP1" ~ "mock community",
    Treatment %in% "PCP2" ~ "mock community",
    Treatment %in% "PCP3" ~ "mock community",
    Treatment %in% "NCP1" ~ "negative control",
    Treatment %in% "NCP2" ~ "negative control",
    Treatment %in% "NCP3" ~ "negative control",
    Treatment %in% "Control" ~ "0 mg/kg metformin",
    Treatment %in% "T2" ~ "25 mg/kg metformin",
    Treatment %in% "T3" ~ "50 mg/kg metformin",
    Treatment %in% "T4" ~ "75 mg/kg metformin"
  ),
  Age = case_when(
    Sample.Date %in% "8/8/19" ~ "40 weeks",
    Sample.Date %in% "10/16/19" ~ "50 weeks",
    Sample.Date %in% "12/20/19" ~ "60 weeks",
    Sample.Date %in% "1/1/2020" ~ "Control"
  )) %>% 
  mutate(TreatmentAge = paste(Treatment, Age, sep = "_"))

## get into SRA submission format
# template: MIMARKS.survey.host-associated.5.0
# https://submit.ncbi.nlm.nih.gov/biosample/template/?package-0=MIMARKS.survey.host-associated.5.0&action=definition

forsra <- data.frame(
  sample_name = sampdf$.sample_name,
  sample_title = paste(sampdf$.sample_name, sampdf$TreatmentAge, sep = "_"),
  organism = "cloacal metagenome",
  collection_date = sampdf$Sample.Date,
  env_broad_scale = "ENVO_00000446",
  env_local_scale = "ENVO_00000002",
  env_medium = "Housed",
  geo_loc_name = "Pennsylvania USA",
  host = "Cobb 500 broiler breeder hen",
  lat_lon = "40.8132368, -77.8698626",
  chem_administration = "metformin supplemented in feed",
  host_age = sampdf$Age,
  host_disease = "none",
  host_sex = "female",
  host_tissue_sampled = "cloacal swab",
  samp_collect_device = "sterile swab",
  description = sampdf$TreatmentAge
  
  
  
)

# write to CSV
write.table(forsra, file = "data/sra_biosample_attributes.txt", sep  = "\t", row.names = FALSE)

### get numbers for males to remove 16S reads
males <- samdat_tbl(ps) %>% 
  anti_join(samdat_tbl(psf))

males$.sample_name
