
## make a phyloseq object with metadata

# Emily Van Syoc 11/2021

# get data
require(dada2)
require(tidyverse)
require(phyloseq)
require(Biostrings)

## get dada2 objects
load("./dada2-data/metforminasv-table.RData")
load("./dada2-data/metformintax-table.RData")

## get sample keys (no metadata right now)
key <- read.table("./data/sampleIDs-with-sequencing-controls.txt", sep = "\t", header = TRUE) %>% 
  column_to_rownames(var = "NovaGeneID")

## get dada2 objects into phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax))

# fix sample IDs
sample_names(ps) <- str_remove(sample_names(ps), ".trim_1P.fastq")

# merge
sample_data(ps) <- key

### from dada2 tutorial: fix ASV names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


# write phyloseq object
save(ps, file = "./data/raw-ps.RData")
