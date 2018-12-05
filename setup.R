

library(SeqArray)
library(tibble)
library(GenomicRanges)
library(tibble)
library(dplyr)
library(pbmcapply)


gff <- "/media/chris/Seagate Expansion Drive/Divergent_Plutella/cleanedGenomes/Aassec_busco/gffs/EOG090W0QXM.gff"
inDir <- "/media/chris/Seagate Expansion Drive/Divergent_Plutella/cleanedGenomes/Aassec_busco/gffs/"
outDir <- "/media/chris/Seagate Expansion Drive/Divergent_Plutella/cleanedGenomes/Aassec_busco/gffs/"

## files for analysis
VCF <- "~/Desktop/Tree-TipR/"
GDS <- "~/Desktop/setupGear/Plutella_filtered_noMAF.recode.gds"
  #"~/Desktop/Tree-TipR/PlutellaSNP.GDS"
df <- readr::read_tsv("/media/chris/Seagate Expansion Drive/whitePupae/GCF_001853355.1_ASM185335v1_genomic.fna.fai", col_names = FALSE)
#get header from VCF will be needed for windows
contigMD <- seqVCF_Header("/media/chris/Seagate Expansion Drive/whitePupae/test2.vcf")$contig
contigMD <- tibble::data_frame(ID = "20", length = 63000000)
contigMD <- tibble::data_frame(ID = df[[1]], length = df[[2]])

contigMD <- filter(contigMD, ID %in% c("NW_017535808.1", "NW_017536897.1", "NW_017537045.1", "NW_017537168.1", "NW_017537234.1"))

# converty VCF to GDS
seqVCF2GDS(vcf.fn = "Dorsalis/pooledsAllSitesDorsalis_filtered.vcf.gz",
             parallel = 6, out.fn = "Dorsalis/pooledsAllSitesDorsalis_filtered.gds",
           fmt.import = c("GT"), storage.option = "ZIP_RA")

# Parameters
windowSize <- 100000
stepSize <- 0
minSites <- 0.25

#open GDS
GDS <- seqOpen("/media/chris/Seagate Expansion Drive/whitePupae/variantsOnly/Allsamples2.gds", allow.duplicate = TRUE)




samples <- SeqArray::seqGetData(gdsfile = GDS, var.name = "sample.id")
pops <- tibble::data_frame(Sample = samples[c(1:6)], Population = c("dorsalis", "dorsalis", "tryoni", "tryoni", "hybrid", "hybrid"))


pops <- tibble::data_frame(Sample = samples[1:100], Population = c(rep("I1", 33), rep("I2", 34), rep("o", 33)))
pops <- tibble::data_frame(Sample = samples[1:6], Population = rep(c("pop1", "pop2"), each = 3))

loci <- windowMaker(contigMD, windowSize,  nCores = 5)

system.time(windowMaker(contigMD, windowSize, stepSize = stepSize, nCores = 5))
## test speed

system.time(geaR::getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 6, pops = pops, ploidy = 2, pairwiseDeletion = TRUE))
# 170MB data
# user  system elapsed
# 0.086   0.062  25.106
test <-geaR::getDiversityStats(GDS, loci[1:1000], minSites = 0.0005, nCores = 5, pops = pops, ploidy = 2, pairwiseDeletion = TRUE)

# Test counting N
test <- getDiversityStats(GDS, loci, minSites = 0.2, nCores = 5, pops = pops, ploidy = 2, countN = FALSE)


dData <- getDiversityStats(GDS, loci, minSites = 0.02, nCores = 4, pops = pops, ploidy = 2, pairwiseDeletion = TRUE)


system.time(getDiversityStats(GDS, loci[1:100], minSites = 0.0005, nCores = 5, pops = pops2, ploidy = 2))



## Codons
fastaName <- "~/Desktop/humaTestData/Homo_sapiens.GRCh38.dna.chromosome.13.fa"
gffName <- "~/Desktop/Pea Project/References/P.xylostella_exon.gff"
gffName <- "~/Desktop/humaTestData/Homo_sapiens.GRCh38.93.chromosome.13.gff3"


genome <- readDNAStringSet(fastaName)
FL <- GRangesList(featureList)
