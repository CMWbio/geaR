

library(SeqArray)
library(tibble)
library(GenomicRanges)
library(tibble)
library(dplyr)
library(pbmcapply)


## files for analysis
VCF <- "~/Desktop/Tree-TipR/"
GDS <- "~/Desktop/setupGear/Plutella_filtered_noMAF.recode.gds"
  #"~/Desktop/Tree-TipR/PlutellaSNP.GDS"

#get header from VCF will be needed for windows
contigMD <- seqVCF_Header("~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz")$contig

# converty VCF to GDS
seqVCF2GDS(vcf.fn = "/media/chris/PhD/Introgression/Plutella_filtered_noMAF.recode.vcf.gz",
           parallel = 6, out.fn = "/media/chris/PhD/Introgression/Plutella_filtered_noMAF.recode.gds",
           storage.option = seqStorageOption(compression = "ZIP_RA.fast", geno.compress = "none"), optimize = TRUE)

## Parameters
windowSize <- 50000
stepSize <- 0
minSites <- 0.005

#open GDS
GDS <- seqOpen("/media/chris/PhD/Introgression/Plutella_filtered_noMAF.recode.gds", allow.duplicate = TRUE)




samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
pops <- data_frame(Sample = samples, Population = c("PxA", "PxA", "Pa", "PxH", "PxH", "PxH",
                                                    "PxH", "PxH", "PxH", "PxH", 'Pa', "PxH",
                                                    "Pa", "PxA", "Pa", "PxA", "Pa", "PxA", "Pa", "Pa",
                                                    "PxA", "PxA", "Pa", "Pa", "PxA", "Pa", "Pa", "Pa", "Pa"))
pops2 <- data_frame(Sample = samples[1:6], Population = rep(c("pop1", "pop2"), each = 10))

loci <- windowMaker(contigMD, windowSize, nCores = 5)

system.time(windowMaker(contigMD, windowSize, nCores = 5))
## test speed

system.time(geaR::getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 6, pops = pops, ploidy = 2, pairwiseDeletion = TRUE))
# 170MB data
# user  system elapsed
# 0.086   0.062  25.106
test <-geaR::getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops, ploidy = 2, pairwiseDeletion = TRUE)

# Test counting N
test <- getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops, ploidy = 2, countN = FALSE)
test <- getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops, ploidy = 2, countN = TRUE)


system.time(getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops2, ploidy = 2))
