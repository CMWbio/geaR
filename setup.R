

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
setwd("/media/chris/Seagate Expansion Drive/whitePupae/")

#"~/Desktop/Tree-TipR/PlutellaSNP.GDS"
df <- readr::read_tsv("GCF_000789215.1_ASM78921v2_genomic.fna.fai", col_names = FALSE)

#get header from VCF will be needed for windows
# contigMD <- tibble::data_frame(ID = df[[1]], length = df[[2]])

# latifrons contigs of interest
contigMD <- dplyr::filter(contigMD, ID %in% c("NW_011876398.1", "NW_011876372.1", "NW_011876235.1"))

# Parametersd
windowSize <- 100000
stepSize <- 0
minSites <- 0.25

loci <- windowMaker(contigMD, windowSize, stepSize = 0, nCores = 5)

which <- GRangesList(GRanges(seqnames = "NW_011876372.1", ranges = IRanges(start = 1, end = 2419413)), GRanges(seqnames = "NW_011876398.1", ranges = IRanges(start = 1, end = 5850278)))
exons <- geaR::getFeatures(gffName = "references/GCF_000789215.1_ASM78921v2_genomic.gff", includeRange = which, feature = )





GDS <- seqOpen("Dorsalis/AllSamplesAllSitesDorsalis_filtered.gds", allow.duplicate = TRUE)

samples <- SeqArray::seqGetData(gdsfile = GDS, var.name = "sample.id")
pops <- tibble::data_frame(Sample = samples[c(1:6)], Population = c("B. dorsalis ", "B. dorsalis ", " B. tryoni ", " B. tryoni ", " hybrid", " hybrid"))

VariantOnlyMandF <- getDiversityStats(GDS, loci, minSites = 0.0002, nCores = 4, pops = pops, ploidy = 2, pairwiseDeletion = TRUE)


# to set up genome
genome <- Biostrings::readDNAStringSet("references/GCF_000789215.1_ASM78921v2_genomic.fna")

## need to change names as whole name is imported
genome@ranges@NAMES <- gsub(" .*$","", genome@ranges@NAMES)

