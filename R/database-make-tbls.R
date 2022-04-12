## add ASV and tax table as database tbl based on project iD lu

# EVS 2/2022

require(tidyverse)
require(phyloseq)
require(microViz)

# get project ID from database
ids <- readxl::read_xlsx(path = "data/database-projectID.xlsx") %>% 
  # this read in as POSIXct
  mutate(Sample_Date = as.Date(Sample_Date) )

# get ASV table
load("data/ps-decontam-filtered-counts.RData")

# get sample data
sampdf <- samdat_tbl(pscount) %>% 
  # convert to POSIXct
  mutate(Sample.Date = str_replace_all(Sample.Date, "19", "2019")) %>%  
  mutate(Sample_Date = as.Date(Sample.Date,
                               tryFormats = c("%m/%d/%Y"))) %>% 
  rename(Animal_ID = Bird) %>% 
  # join to ids to get Project-ID
  full_join(ids) %>% 
  column_to_rownames(var = ".sample_name")


# add to phyloseq 
sample_data(pscount) <- sampdf

# change sample names to Project ID
sample_names(pscount) <- sampdf$Project_ID

## ---- clean and make vertical for database

# get ASV table
asv <- otu_get(pscount) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Project_ID") %>% 
  # make vertical
  pivot_longer(cols = starts_with("ASV"), names_to = "ASV", values_to = "count")
  

# write
write.table(asv, file = "data/for-database/tbl-ChickenMetformin_ASVs.txt", sep = "\t", row.names = FALSE)

# get tax table
tax <- tax_table(pscount) %>% 
  as.data.frame() %>% 
  # make vertical
  rownames_to_column(var = "ASV")

# write.table
write.table(tax, file = "data/for-database/tbl-ChickenMetformin_tax_table.txt",
           sep = "\t", row.names = FALSE)
