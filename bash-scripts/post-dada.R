## finishing the dada2 pipeline
# post -dada

## --- CHANGE THESE VARIABLES ----

# set working directory
OUTPATH <- "/gpfs/group/evk5387/default/Novogene/metformin"

# load previous RData
load(paste0(OUTPATH, "/all-dada2-objects.RData"))

## path to reference database
# downloaded from: https://benjjneb.github.io/dada2/training.html
DB = "/gpfs/group/evk5387/default/databases/silva_nr99_v138.1_train_set.fa.gz"

### --- CODE ----

require(dada2)

# merge paired reads
mergers <- mergePairs(dadaF = dadaForward,
                      derepF = forward,
                      dadaR = dadaReverse,
                      derepR = reverse,
                      verbose = TRUE)
					  
# save intermediate progress
save.image(paste0(OUTPATH, "/working-image.RData"))

# construct sequence table of ASVs
seqtab <- makeSequenceTable(samples = mergers)

# save intermediate progress
save.image(paste0(OUTPATH, "/working-image.RData"))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(unqs = seqtab, 
                                    method = "consensus",
                                    verbose = TRUE)
									
# save intermediate progress
save.image(paste0(OUTPATH, "/working-image.RData"))

# assign taxonomy using the Silva database
tax <- assignTaxonomy(seqs = seqtab.nochim, 
                      refFasta = DB, 
                      verbose = TRUE)
					  
# save intermediate progress
save.image(paste0(OUTPATH, "/working-image.RData"))

## This is the end of the dada2 algorithm

## ---- WRITE OUTPUT ----

# write ASV table to file
write.table(seqtab.nochim, file = paste0(OUTPATH, "asv-table.txt"), sep = "\t", row.names = FALSE)

# save ASB table as RData object
save(seqtab.nochim, file = paste0(OUTPATH, "asv-table.RData"))

# write taxonomy table to file
write.table(tax, file = paste0(OUTPATH, "tax-table.txt"), sep = "\t", row.names = FALSE)

# save tax table as RData object
save(tax, file = paste0(OUTPATH, "tax-table.RData"))

# save all objects
save.image(file = paste0(OUTPATH, "/all-dada2-objects.RData"))

## ---- track reads through the pipeline ----

# define function (from dada2 tutorial)
getN <- function(x) sum(getUniques(x))

# create dataframe
track <- cbind(sapply(dadaForward, getN), sapply(dadaReverse, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# change names
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- fwdNames

# write to file
write.table(track, file = paste0(OUTPATH, "track-reads-dada.txt"), sep = "\t") 

# print finished message
cat("done!")


