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
#' @param ploidy \code{numeric} number of alleles
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


getDiversityStats <- function(GDS, loci, minSites = 0.5, nCores = 1, pops = NULL, stats = "all", ploidy = 2){
  #browser()
  if(length(pops) == 0){
    samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
    pops <- data_frame(Sample = samples, Population = samples)
  }

  popList <- split(pops, pops$Population)

  nam <- names(popList)
  pairs <- outer(nam, nam, paste, sep = "///")
  pairs <- lapply(pairs, function(z){
    pair <- unlist(strsplit(z, split = "///"))
    if(pair[1] == pair[2]) pair <- NULL
    pair <- sort(pair)
    return(pair)
  })
  pairs <- pairs[!duplicated(pairs)]
  pairs <- Filter(Negate(is.null), pairs)


  stopifnot(minSites < 1 & minSites > 0)
  ## Iterating through the window list by scaffold
  div <- pbmclapply(seq(length(loci)), mc.cores = 5, function(locusN){
    locus <- loci[[locusN]]

    genoMat <- getGenotypes(GDS = GDS, locus = locus, minSites = minSites, nucleotide = FALSE, ploidy = ploidy, pops = pops)

    distMat <- genoDist(genoMat)

    seqname <- locus@seqnames@values
    start <- locus@ranges@start
    end <- start + locus@ranges@width
    windowMid <- (start + end) /2
    snpMid <- floor(mean(as.numeric(rownames(genoMat)), na.rm = TRUE))
    nSites <- ifelse(length(genoMat), nrow(genoMat), 0)

    dxy <-  neisDxy(distMat, popList, pairs, ploidy = 2)
    pi <- neisPi(distMat, popList, ploidy = 2)
    da <- neisDa(dxy, pi)



    bind_cols(data_frame(SeqName = seqname, Start = start, End = end, windowMid, snpMid, nSites), pi, dxy, da)

  })

  div <- bind_rows(div)

  return(div)
}

