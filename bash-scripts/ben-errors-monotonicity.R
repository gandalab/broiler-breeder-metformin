## R script to run dada2 algorithm on trimmed paired-end reads
## with monotonicity on binned quality scores from Ben Callahan

## ---- SET VARIABLES ----

## modify these variables:

# path to filtered and cleaned reads
CLEANEDPATH = "/gpfs/group/evk5387/default/Novogene/metformin/trimmed"

# path to output tables
OUTPATH = "/gpfs/group/evk5387/default/Novogene/metformin"

# database to silva training set
# downloaded from: https://benjjneb.github.io/dada2/training.html
DB = "/gpfs/group/evk5387/default/databases/silva_nr99_v138.1_train_set.fa.gz"

# paired end patterns 
FILTEREDF = ".trim_1P.fastq"
FILTEREDR = ".trim_2P.fastq"

# ---- install and read data ----

### packages must be previously installed
require(dada2)
require(tidyr)
require(phyloseq)
require(magrittr)

## test that pathway works
#if(!list.files(CLEANEDPATH)) {
#  cat("Can't read file pathway or files are not present")
#}

## ---- core dada algorithm ----

# get forward and reverse reads
forward <- sort(list.files(CLEANEDPATH, pattern = FILTEREDF, full.names = TRUE))
reverse <- sort(list.files(CLEANEDPATH, pattern = FILTEREDR, full.names = TRUE))

# check to make sure that the lengths of both files are the same and that they match
fwdNames <- sapply(strsplit(basename(forward), FILTEREDF), `[`, 1)
revNames <- sapply(strsplit(basename(reverse), FILTEREDR), `[`, 1)

# error catch
if(length(fwdNames) != length(revNames)) {
  stop("The number of forward and reverse files do not match.")
} else {
  
  if(any(!fwdNames%in% revNames)) {
    
    stop("Forward and reverse reads are out of order.")
  }
}


# perform error learning
errF <- learnErrors(forward, 
                    multithread = TRUE,
					verbose = TRUE)
errR <- learnErrors(reverse, 
                    multithread = TRUE,
					verbose = TRUE)
					
# enforce monotonicity
make.monotone.decreasing <- function(v) sapply(seq_along(v), function(i) max(v[i:length(v)]))
errF.md <- t(apply(getErrors(errF), 1, make.monotone.decreasing))
errR.md <- t(apply(getErrors(errR), 1, make.monotone.decreasing))

# save progress
save.image(file = paste0(OUTPATH, "/monotonicity-error.RData"))

# perform denoising on forward and reverse reads
dadaForward <- dada(derep = forward, 
                    err = errF.md, 
                    multithread = TRUE,
					verbose = TRUE)
dadaReverse <- dada(derep = reverse, 
                    err = errR.md, 
                    multithread = TRUE,
					verbose = TRUE)

## save error plots in PDF format
# hacky fix to plot errors
# forward
errF.md.full <- errF
errF.md.full$err_out <- errF.md
dimnames(errF.md.full$err_out) <- dimnames(as.matrix(errF$err_out))
# reverse
errR.md.full <- errR
errR.md.full$err_out <- errR.md
dimnames(errR.md.full$err_out) <- dimnames(as.matrix(errR$err_out))


# plot forward errors
pdf(paste0(OUTPATH, "/forward-errorplot.pdf"))
plotErrors(errF.md.full, nominalQ = TRUE)
dev.off()

# plot reverse errors
pdf(paste0(OUTPATH, "/reverse-errorplot.pdf"))
plotErrors(errR.md.full, nominalQ = TRUE)
dev.off()

# save progress
save.image(file = paste0(OUTPATH, "/all-dada2-objects.RData"))

