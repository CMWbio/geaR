

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
df <- readr::read_tsv("/media/chris/Seagate Expansion Drive/whitePupae/GCF_000789215.1_ASM78921v2_genomic.fna.fai", col_names = FALSE)
#get header from VCF will be needed for windows
contigMD <- tibble::data_frame(ID = "20", length = 63000000)
contigMD <- tibble::data_frame(ID = df[[1]], length = df[[2]])

contigMD <- filter(contigMD, ID %in% c("NW_011876398.1", "NW_017536897.1", "NW_017537045.1", "NW_017537168.1", "NW_017537234.1"))

# converty VCF to GDS
seqVCF2GDS(vcf.fn = "Dorsalis/WholeGenomeAnalysis/QTLAllSitesDorsalis_withBo_filtered.vcf.gz",
             parallel = 6, out.fn = "Dorsalis/QTLAllSitesDorsalis_withBo_filtered.gds",
           fmt.import = c("GT"), storage.option = "ZIP_RA")

# Parameters
windowSize <- 100000
stepSize <- 0
minSites <- 0.25

#open GDS
setwd("/media/chris/Seagate Expansion Drive/1000Genomes/")



#get header from VCF will be needed for windows
contigMD <- tibble::data_frame(ID = "chr1", length = 10000000)

# Parametersd
windowSize <- 10000
stepSize <- 0
minSites <- 0.25

loci <- windowMaker(contigMD, windowSize, stepSize = 0, nCores = 5)

which <- GRangesList(GRanges(seqnames = "NW_011876372.1", ranges = IRanges(start = 1, end = 650000)), GRanges(seqnames = "NW_011876398.1", ranges = IRanges(start = 3800000, end = 5850278)))
exons <- geaR::getFeatures(gffName = "/media/chris/Seagate Expansion Drive/whitePupae/references/GCF_000789215.1_ASM78921v2_genomic.gff", nCores = 4, includeRange = which, feature = "gene:cds")





GDS <- seqOpen("ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.gds", allow.duplicate = TRUE)

samples <- SeqArray::seqGetData(gdsfile = GDS, var.name = "sample.id")
pops <- tibble::data_frame(Sample = samples[c(1:1000)], Population = rep(c("P1", "P2", "P3", "P4"), each = 250))
pops2 <- split(pops, pops$Population)



VariantOnlyMandF <- getDiversityStats(GDS, loci, minSites = 0.01, nCores = 4, stats = "dxy", pops = pops, ploidy = 2, pairwiseDeletion = TRUE)


# to set up genome
genome <- Biostrings::readDNAStringSet("references/GCF_000789215.1_ASM78921v2_genomic.fna")

## need to change names as whole name is imported
genome@ranges@NAMES <- gsub(" .*$","", genome@ranges@NAMES)

