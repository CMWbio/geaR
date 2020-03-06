# Introduction

The package geaR is designed to aid with evolutionary and population genetic analysis across whole genome genotype data in the GDS format.
Analysis can be carried out using classical approaches such as sliding windows, however the main strength of geaR is the ability to include only ceratin types of features in the analysis.

This package is currently under development, use at own risk. 

# Installation

```

install.packages("BiocManager")
install.packages("devtools")
devtools::install_github("CMWbio/geaR")

```
# Quick usage 
For detailed usage please see here.
The easiest way to use geaR is to use the cog/gear object interface. 
This example will use a reasonably sized human chromosome 20 downloaded from here[ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/] on a 4 core machine.


```
## Make 100kb windows across the genome.
### This requires a dataframe of scaffold/chromosome lengths.

library(geaR)
library(tibble)
library(SeqArray)

### if you dont know the lengths, it can easily be determined from your reference index (.fai)
chr20_df <- tibble(ID = "chr20", length = 64444167)
loci <- windowMaker(chr20_df, windowSize = 100000, stepSize = 0,  nCores = 4)

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




```
