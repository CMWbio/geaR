#' Calculate diversity statistics
#'
#' @description Imports genotypes and calculates diversity formula
#'
#' @details Authour: Chris Ward
#' Uses a \code{list} containing \code{GRanges} or a \code{GRangesList} to import genotypes. Diversity statistics are then calculated.
#'
#' @param GDS \code{GDS} object with variant data to import genotypes from.
#' @param loci \code{list} or \code{GRangesList}. Loci to import genotypes for.
#' @param minSites \code{numeric} minimum number of sites as a proportion of loci length. Default 0.5 (ie 50%)
#' @param nCores \code{numeric} number of cores to run in parallel
#' @param pops \code{data_frame} containing Sample ID and poplations. Default is to treat individuals as populations \code{c("none")}
#' @param stats \code{character}. Vector containing all diversity stats to calculate. default \code{c("all")}
#'
#'
#' @return A \code{data_frame} of selected Diversity statistics
#'
#'
#' @examples
#'
#'
#'
#' @export
#' @rdname getDiversityStats


getDiversityStats <- function(GDS, loci, minSites = 0.5, nCores = 1, pops = "none", stats = "all", ploidy = 2){

  stopifnot(minSites < 1 & minSites > 0)
   popList <- split(pops, pops$Population)
    ## Iterating through the window list by scaffold
    system.time(div <- pbmclapply(seq(length(loci)), mc.cores = 5, function(locus){

      genoMat <- getGenotypes(GDS = GDS, locus = loci[[locus]], minSites = 0.00005, nucleotide = FALSE)

      distMat <- genoDist(genoMat)

      if(pops == "none"){
        samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
        pops <- data_frame(Sample = samples, Population = samples)
      }

     neisDxy(distMat, popList, ploidy)

   }))




}

