

library(SeqArray)
library(tibble)
library(GenomicRanges)
library(tibble)
library(dplyr)
library(pbmcapply)

######## new version with S4 methods
arg <- makeCog(analysisType = "args", ploidy = 2, nCores = 4, minSites = 0.02, pairwiseDeletion = TRUE, removeIndels = TRUE)
diversity <- makeCog(analysisType = "diversityFULL", stats = "FTD")
admix <- makeCog(analysisType = "admixture", threePop = list(c("tryoni", "hybrid", "dorsalis")),
                     fourPop = rep(list(c("hybrid", "tryoni", "dorsalis", "oleae")), 2))
outloci <-  makeCog(analysisType = "outputLoci", outputDirectory = "~/Desktop/setupGear/", alleles = "seperate", removeIndels = TRUE)
outTrees <-  makeCog(analysisType = "outputTrees", outputDirectory = "~/Desktop/setupGear/", alleles = "seperate", removeIndels = TRUE)


gear <- makeGear(loci, populations = pops, outgroup = "oleae", cogs = list(diversity, arg))

gear <- analyzeGear(GDS = AllSamplesVariant, gear)

######
CDS <- makeFeatures("/media/chris/Seagate Expansion Drive/whitePupae/references/GCF_000347755.3_Ccap_2.1_genomic.gff", feature = "gene:cds", nCores = 4, longestIsoform = TRUE)
exon <- makeFeatures("/media/chris/Seagate Expansion Drive/whitePupae/references/GCF_000347755.3_Ccap_2.1_genomic.gff", feature = "gene:exons", nCores = 4)



##### testing read.me code

library(geaR)
library(tibble)
library(SeqArray)

### if you dont know the lengths, it can easily be determined from your reference index (.fai)
chr20_df <- tibble(ID = "chr20", length = 64444167)
loci <- windowMaker(chr20_df, windowSize = 1000000, stepSize = 0,  nCores = 4)

## Next we will construct our analysis
### first loading the GDS and setting up population definitions.
GDS <- seqOpen("/media/chris/Seagate Expansion Drive/geaR_Test/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.gds")

### We can easily query the GDS for sample names.
samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
### for the example populations are arbitrary so we will just set 100 equal populations each with 10 samples
pops <- tibble(Sample = samples[1:1000], Population = rep(paste0("P", 1:100), each = 10))

### general arguments for the analysis
argCog <- makeCog(analysisType = "args", ploidy = 2, nCores = 4, minSites = 0.002, pairwiseDeletion = TRUE, removeIndels = TRUE)

### set up diverstity analysis
divCog <- makeCog(analysisType = "diversityFULL", stats = "all")

### set up admixture cog
admixCog <- makeCog(analysisType = "admixture", threePop = list(),
                 fourPop = list("all"))

### set up cog to output loci to file

outlociCog <-  makeCog(analysisType = "outputLoci", outputDirectory = "~/Desktop/setupGear/", alleles = "seperate", removeIndels = TRUE)


### build the gear object for analysis
### We will also arbitrarily define P100 as the outgroup for this analysis
gear <- makeGear(loci[1:20], populations = pops, outgroup = "P100", cogs = list(argCog, divCog, outlociCog))

### Run the three analyses 
gear <- analyzeGear(GDS, gear)


############ codons 

GDS <- seqOpen("/media/chris/Seagate Expansion Drive/whitePupae/Dorsalis/AllSamplesAllSitesDorsalis_filtered.gds", allow.duplicate = TRUE)
samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
### for the example populations are arbitrary so we will just set 100 equal populations each with 10 samples
pops <- tibble::tibble(Sample = samples, Population = samples)
cds <- makeFeatures("/media/chris/Seagate Expansion Drive/whitePupae/references/GCF_000789215.1_ASM78921v2_genomic.gff", feature = "gene:cds", nCores = 7, geneIDField = "Name", longestIsoform = TRUE)

genome <- Biostrings::readDNAStringSet("/media/chris/Seagate Expansion Drive/whitePupae/references/GCF_000789215.1_ASM78921v2_genomic.fna")

## need to change names as whole name is imported
genome@ranges@NAMES <- gsub(" .*$","", genome@ranges@NAMES)

tic()
fourF <- buildCodonDB(genome, exons = cds[1:10], sqlDir = NULL)
toc()
tic()
f <- validate4FoldCodons(GDS, input = fourF, nCores = 4, pops = pops)
toc()

