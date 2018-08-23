

library(SeqArray)
library(tibble)
library(GenomicRanges)


## files for analysis
VCF <- "~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz"
GDS <- "~/Desktop/Tree-TipR/PlutellaSNP.GDS"

#get header from VCF will be needed for windows
contigMD <- seqVCF_Header("~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz")$contig

# converty VCF to GDS
seqVCF2GDS(vcf.fn = "~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz",
           parallel = 6, out.fn = "~/Desktop/Tree-TipR/PlutellaSNP.GDS")

## Parameters
windowSize <- 100000
stepSize <- 100000
minSites <- 1000

#open GDS
gds <- seqOpen(GDS)

genotypes <-  seqGetData(gdsfile = gds, var.name = "genotype")

