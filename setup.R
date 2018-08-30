

library(SeqArray)
library(tibble)
library(GenomicRanges)
library(tibble)
library(dplyr)
library(pbmcapply)


## files for analysis
VCF <- "~/Desktop/Tree-TipR/"
GDS <- "~/Desktop/Tree-TipR/PlutellaSNP.GDS"

#get header from VCF will be needed for windows
contigMD <- seqVCF_Header("~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz")$contig

# converty VCF to GDS
seqVCF2GDS(vcf.fn = "~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz",
           parallel = 6, out.fn = "~/Desktop/Tree-TipR/PlutellaSNP.GDS",storage.option = "ZIP_RA", optimize = TRUE)

## Parameters
windowSize <- 100000
stepSize <- 0
minSites <- 0.005

#open GDS
GDS <- seqOpen(GDS, allow.duplicate = TRUE)


pops <- data_frame(Sample = samples[1:6], Population = rep(c("pop1", "pop2"), each = 3))
pops2 <- data_frame(Sample = samples[1:20], Population = rep(c("pop1", "pop2"), each = 10))

samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")

loci <- windowMaker(contigMD, windowSize, nCores = 5)

system.time(windowMaker(contigMD, windowSize, nCores = 5))
## test speed

system.time(getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops, ploidy = 2))
# user  system elapsed
# 0.084   0.031   6.691
test <- getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops, ploidy = 2)

# Test counting N
test <- getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops, ploidy = 2, countN = FALSE)
test <- getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops, ploidy = 2, countN = TRUE)


system.time(getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops2, ploidy = 2))
