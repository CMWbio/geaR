###################################################################################################
## Libraries
###################################################################################################
library(SeqArray)
library(magrittr)
library(tidyverse)
library(pbmcapply)
library(ggplot2)

###################################################################################################
## Test data
###################################################################################################
## File paths to objects
VCF <- "~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz"
GDS <- "~/Desktop/Tree-TipR/PlutellaSNP.GDS"

## 1000 Genomes VCF - chr1-2
seqVCF2GDS(vcf.fn = VCF, out.fn = VCF, parallel = 6, verbose = TRUE)
genofile <- seqOpen(gds.fn = GDS)
samp.id <- seqGetData(genofile, "sample.id")
seqSetFilter(object = genofile, sample.id=samp.id[c(1:10)])

## Mock population structure
pops <- data.frame(samples = seqGetData(genofile, "sample.id"),
                   populations = c(rep("pop1",2),rep("pop2",2),rep("pop3",2),rep("pop4",2),rep("pop5",2)))

## Parameters
windowSize <- 100000
step <- 100000
minSites <- 1000

###################################################################################################
## Generating windows
###################################################################################################
getWindows <- function(path_vcf, windowSize, step, chromosomes) {

  # browser()

  ## Building chromosome dataframe
  vcfHeader <- VariantAnnotation::scanVcfHeader(file = path_vcf)@header$contig
  vcfHeader <- as.data.frame(vcfHeader)
  vcfHeader$chr <- rownames(vcfHeader)
  vcfHeader
  vcfHeader <- filter(.data = vcfHeader, chr == chromosomes)

  ## Building GRanges object of chromosomes/scaffolds
  contig_granges <- select(.data = vcfHeader, chr, end = length)
  contig_granges <- mutate(.data = contig_granges, start = 1, strand = ".", end = end)
  contig_granges <- select(.data = contig_granges, chr,start, end, strand)
  contig_granges <- GenomicRanges::makeGRangesFromDataFrame(contig_granges)

  ## List of windows
  window_list <- lapply(chromosomes, function(x){

    ## Subsetting GRange object by chromosome
    window_chr <- contig_granges[GenomicRanges::seqnames(contig_granges) == x]

    ## Splitting into windows
    windows_chr <- GenomicRanges::slidingWindows(x = window_chr, width = windowSize, step = step)
    windows_chr <- unlist(windows_chr) ## Need to unlist the slidingWindow object

    ## Turning window-GRanges into GRangesList
    windows_chr <- as(windows_chr,"GRangesList")
    windows_chr

  })

  ## Returning list object
  names(window_list) <- chromosomes
  window_list

}

## Creating window list
system.time({
  window_list <- getWindows(path_vcf = VCF, chromosomes = 1,windowSize = windowSize, step = windowSize)
})

###################################################################################################
## Genotype Matrices
###################################################################################################
getGenotypeMatrices <- function(gds_file, windows_list, minSites = minSites, nCores){

  ## Iterating through the window list by scaffold
  by_scaffold <- lapply(names(windows_list), function(scaffoldName){

    ## Iterating through windows for current scaffold
    by_window <- pbmclapply(windows_list[[scaffoldName]], mc.cores = nCores, function(currentWindow){

      ## Subsetting GDS file
      seqSetFilter(object = gds_file, currentWindow) ## Setting filter (i.e. window to get variants from)
      arr <- seqGetData(gdsfile = gds_file, var.name = "genotype")

      ## Filtering empty arrays/arrays with too few variants
      if(dim(arr)[3] < minSites || is.null(arr)) {

        print("Incompatible window - not enough information")
        return(NULL)

      } else {

        ## Number of variants in window
        varNumber <- length(seqAlleleCount(gds_file))

        ## Building window identifier
        chr <- as.vector(GenomicRanges::seqnames(currentWindow)@values) ## Getting the current scaffold
        # chr <- paste0("chr",chr)
        start <- GenomicRanges::start(currentWindow)
        end <- GenomicRanges::end(currentWindow)
        win_id <- paste(chr,start,end,sep = "_")

        ## Getting sample and variant IDs for col/row nameing
        sample.id <- seqGetData(gdsfile = gds_file, var.name = "sample.id")
        variant.id <- seqGetData(gdsfile = gds_file, var.name = "variant.id")

        ## Duplicating row-names (duce for splitting into haplotypes)
        sample.id.1 <- paste0(sample.id,"_1")
        sample.id.2 <- paste0(sample.id,"_2")

        ## Combining sample ID with window ID information
        s_id_1 <- paste(sample.id.1,win_id,varNumber, sep = "_")
        s_id_2 <- paste(sample.id.2,win_id,varNumber, sep = "_")

        ## Row order to order matrix by
        sample.order <- c(rbind(s_id_1,s_id_2)) ## new naming convention - same order
        # sample.order <- c(rbind(sample.id.1,sample.id.2)) ## Old naming convention - no window information

        ## Building matrix 1
        mat.1 <- matrix(arr[1,,],
                        nrow=dim(arr)[2],
                        ncol=dim(arr)[3], dimnames = list(s_id_1,variant.id)
        )
        mat.1[mat.1 == "NA/NA"] <- NA

        ## Building matrix 2
        mat.2 <- matrix(arr[2,,],
                        nrow=dim(arr)[2],
                        ncol=dim(arr)[3], dimnames = list(s_id_2,variant.id)
        )
        mat.2[mat.2 == "NA/NA"] <- NA

        # ## Joining matrices
        mat <- rbind(mat.1,mat.2)

        ## Ordering by rownames
        mat[order(match(rownames(mat),sample.order)),]

      }

    })

  })

  ## Assigning scaffold names to list elements
  names(by_scaffold) <- names(windows_list)

  ## Returning named list of genotype matrices (split into haplotypes) by scaffold
  by_scaffold

}

## Creating genotype matrices list
system.time({
  genotype_matrices <- getGenotypeMatrices(gds_file = genofile, windows_list = window_list, minSites = minSites, nCores = 6)
})

###################################################################################################
## Dxy
###################################################################################################
dxy_neiDist <- function(genotypeMatrices,population_df){

  lapply(genotypeMatrices, function(scaffold){

    by_scaff <- lapply(scaffold, function(mat){

      if(!is.null(mat)){

        ## Genetic distance calculation on haplotype matrix
        dist <- poppr::nei.dist(mat)

        ## Temp for snake data
        # dist <- poppr::provesti.dist(mat)

        ## Converting dist object to matrix
        dist_mat <- as.matrix(dist)
        dist_mat

        ## Window value
        window <- gsub("^.+_.+_(.*_.*_.*_.*)$","\\1",colnames(dist_mat)[1])
        ## NOTE: This will result in out of range error down below
        ## if scaff/chromosome naming is not just numbers (i.e. 1 == good, chr1 == bad)
        ## Quick solution, need to develop

        ## Updating population dataframe
        nm1 <- paste(population_df$samples,"1",window,sep = "_")
        nm2 <- paste(population_df$samples,"2",window,sep = "_")
        nm <- c(rbind(nm1,nm2)) ## Naming populations by their window
        pop <- slice(.data = population_df, rep(1:n(), each = 2)) ## Getting right number of populations for new df
        pop <- as.vector(as.character(pop$pop))

        pop <- data.frame(samples = nm, pops = pop)
        pop <- split(pop, pop$pops)
        dxy <- lapply(pop, function(pop1){

          pop_by_pop <- lapply(pop, function(pop2){
            if(!setequal(pop1,pop2)){
              ##Subsetting distance matrix by
              dist_mat <- dist_mat[as.vector(as.character(pop1$samples)), as.vector(as.character(pop2$samples))]
              df <- data_frame(mean(dist_mat),sd(dist_mat))
              colnames(df) <- paste0(pop1$pops[1], "/" , pop2$pops[1], c("_dxy", "_SDdxy"))
              df
            }
          })

          ## Inidivdual population list elements combined into one row
          bind_cols(pop_by_pop)
        })

        ## All popuilations by window collapsed into one row representing that window with all population comparisons
        dxy <- bind_cols(dxy)
        dxy$window <- window
        dxy <- select(dxy, window, everything())
        dxy

      }

    })

    ## binding windows by scaffold
    bind_rows(by_scaff)

  })

}

## Dxy values by chr/scaffold
system.time({
  dxy <- dxy_neiDist(genotypeMatrices = genotype_matrices, population_df = pops)
})
