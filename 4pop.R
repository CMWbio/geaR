#' 4 pop outgroup statistics
#'
#' @description Imports genotypes across a locus into matrix
#'
#' @details Authours: Chris Ward & Alastair Ludington
#' Uses a GRanges locus to import genotypes, either nucleotide or RAW, from a GDS file
#'
#'
#' @param GDS \code{GDS} object with variant data to import genotypes from
#' @param locus \code{GRanges} Locus to import genotypes for
#' @param minSites \code{numeric} minimum number of sites as a proportion of loci length
#' @param nucleotide \code{logical} Import RAW genotypes or nucleotides
#' @param ploidy \code{numeric} ploidy of sample
#' @param pops \code{data_frame} populaiton dataFrame
#' @param removeIndels removes indels
#'
#'
#' @return A \code{matrix} of genotypes
#'
#'
#' @import SeqArray
#' @import tidyr
#' @import dplyr
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges width
#' @examples
#'
#'
#'
#' @rdname getGenotypes


Tpop <- function(x){




GDS <- seqOpen("Dorsalis/WholeGenomeAnalysis/AllSitesDorsalis_biallelic.vcf.gz", allow.duplicate = TRUE)


if(length(pops) == 0){
  samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
  pops <- data_frame(Sample = samples, Population = samples)
}

popList <- split(pops, pops$Population)

stopifnot(minSites < 1 & minSites > 0)
## Iterating through the window list by scaffold



div <- mclapply(seq(length(loci)), mc.cores = 4, function(locusN){
  locus <- loci[[locusN]]


 AF_list <- lapply(popList, function(x){

      samples <- x$Sample
      seqSetFilter(object = GDS, sample.id = samples)



    if(length(locus)){
      seqSetFilter(object = GDS, variant.sel = locus)
    }


   seqAlleleFreq(GDS)


  })

 AF_df <- bind_cols(AF_list)
 AF_df <- AF_df[complete.cases(AF_df),]

 popString <- c(A = "tryoni", B = "hybrid", C = "dorsalis", D = "oleae")



 f4 <- apply(AF_df, 1, function(y){
   p1 <- y[popString["A"]]
   p2 <- y[popString["B"]]
   p3 <- y[popString["C"]]
   p4 <- y[popString["D"]]

   q1 <- 1 - p1
   q2 <- 1 - p2
   q3 <- 1 - p3
   q4 <- 1 - p4


   # f4 <- (p1 - p2) * (p3 - p4)
   #
   # cor <- (q1 - q2) * (q3 - q4)

   # from simon martin script same result
   f4 <- q1*p2*p3*q4 - p1*q2*p3*q4
   cor <- p1*q2*q3*p4 - q1*p2*q3*p4


   f4 + cor
 })

 denom <- apply(AF_df, 1, function(y){



   p1 <- y[popString["A"]]
   p2 <- y[popString["B"]]
   p3 <- y[popString["C"]]
   p4 <- y[popString["D"]]

   q1 <- 1 - p1
   q2 <- 1 - p2
   q3 <- 1 - p3
   q4 <- 1 - p4

   # multiply by logical to negate lower allele for pd
   pd <- p2* (p2>p3) + p3*(p3>=p2)
   qd <- 1 - pd

   fd <- q1*pd*pd*q4 - p1*qd*pd*q4
   cor <- p1*qd*qd*p4 - q1*pd*qd*p4

   fd <- fd + cor

 })

fd <- sum(f4) / sum(denom)

seqname <- locus@seqnames@values[1]
start <- locus@ranges@start
end <- start + locus@ranges@width

start <- min(start)
end <- max(end)
windowMid <- (start + end) /2

position <- seqGetData(gdsfile = GDS, var.name = "position")
snpMid <- floor(median(position))
nSites <- length(position)

data_frame(SeqName = seqname, Start = start, End = end, windowMid, snpMid, nSites, fd)

})


fd <- bind_rows(div)


library(readr)
write_tsv(fd, "Figure for paper/fd_100kb.tsv")





}



