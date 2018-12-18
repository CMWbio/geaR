#' Calculate diversity statistics
#'
#' @description Imports genotypes and calculates diversity formula
#'
#' @details Authour: Chris Ward
#' Uses a \code{list} containing \code{GRanges} or a \code{GRangesList} to import genotypes. Diversity statistics are then calculated.
#'
#' @param GDS \code{GDS} object with variant data to import genotypes from.
#' @param loci \code{list} or \code{GRangesList}. Loci to import genotypes for.
#' @param minSites \code{numeric} minimum number of sites as a proportion of loci length. Default 0.5 ie 50 percent
#' @param nCores \code{numeric} number of cores to run in parallel
#' @param pops \code{data_frame} containing Sample ID and poplations. Default is to treat individuals as populations \code{c("none")}
#' @param stats \code{character}. Vector containing all diversity stats to calculate. default \code{c("all")}
#' @param ploidy \code{numeric} number of alleles
#' @param pairwiseDeletion \code{numeric} nshould Ns be removed from pairwise comparisons
#'
#'
#' @return A \code{data_frame} of selected Diversity statistics
#'
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#'
#' @export
#' @rdname getDiversityStats


getDiversityStats <- function(GDS, loci, minSites = 0.5, nCores = 1, pops = NULL, stats = "all", ploidy = 2, pairwiseDeletion){
  #browser()



  if(length(pops) == 0){
    samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
    pops <- data_frame(Sample = samples, Population = samples)
  }

  popList <- split(pops, pops$Population)

  nam <- names(popList)
  pairs <- outer(nam, nam, paste, sep = "///")
  pairs <- lapply(pairs,function(z){
    pair <- unlist(strsplit(z, split = "///"))
    if(pair[1] == pair[2]) pair <- NULL
    pair <- sort(pair)
    return(pair)
  })
  pairs <- pairs[!duplicated(pairs)]
  pairs <- Filter(Negate(is.null), pairs)


  stopifnot(minSites < 1 & minSites > 0)
  ## Iterating through the window list by scaffold
  div <- pbmclapply(seq(length(loci)), mc.cores = nCores, function(locusN){
    locus <- loci[[locusN]]

    genoMat <- getGenotypes(GDS = GDS, locus = locus, minSites = minSites, nucleotide = FALSE, ploidy = ploidy, pops = pops)

    distMat <- genoDist(genoMat, pairwiseDeletion)

    seqname <- locus@seqnames@values[1]
    start <- locus@ranges@start
    end <- start + locus@ranges@width

    start <- min(start)
    end <- max(end)
    windowMid <- (start + end) /2

    snpMid <- floor(mean(as.numeric(rownames(genoMat)), na.rm = TRUE))
    nSites <- ifelse(length(genoMat), nrow(genoMat), 0)

    if(length(distMat)){
      dxy <- neisDxy(distMat, popList, pairs, ploidy = ploidy)
      pi <- neisPi(distMat, popList, ploidy = ploidy)
      da <- neisDa(dxy, pi)
      fst <- weightedFst(distMat, popList, pairs, ploidy = ploidy)
    } else{
      dxy <- c()
      pi <- c()
      da <- c()
      fst <- c()
    }


    if("gene" %in% colnames(locus@elementMetadata)) {
      gNames <- paste(unique(locus$gene), collapse = ",")
      bind_cols(data_frame(SeqName = seqname, Start = start, End = end, Gene = gNames, windowMid, snpMid, nSites), pi, dxy, da, fst)}
    else bind_cols(data_frame(SeqName = seqname, Start = start, End = end, windowMid, snpMid, nSites), pi, dxy, da, fst)



  })

  dplyr::bind_rows(div)

}

