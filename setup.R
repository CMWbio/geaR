

library(SeqArray)
library(tibble)
library(GenomicRanges)


## files for analysis
VCF <- "~/Desktop/Tree-TipR/"
GDS <- "~/Desktop/Tree-TipR/PlutellaSNP.GDS"

#get header from VCF will be needed for windows
contigMD <- seqVCF_Header("~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz")$contig

# converty VCF to GDS
seqVCF2GDS(vcf.fn = "~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz",
           parallel = 6, out.fn = "~/Desktop/Tree-TipR/PlutellaSNP.GDS")

## Parameters
windowSize <- 100000
stepSize <- 100000
minSites <- 0.005

#open GDS
GDS <- seqOpen(GDS, allow.duplicate = TRUE)

genotypes <-  seqGetData(gdsfile = gds, var.name = "genotype")

pops <- data_frame(Sample = samples[1:20], Population = rep(c("pop1", "pop2"), each = 10))
samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")


loci <- windowMaker(contigMD, windowSize, nCores = 5)
