#' Calculates maximum genetic distance
#'
#' @description Calculates maximum population genetic distance from hamming distance
#'
#' @details 
#' Calculates maximum population genetic distance from hamming distance
#'
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param pairs list of populaiton pairs generated in \code{getDiversityStats}
#' @param ploidy \code{numeric} number of chromosomes
#'
#'
#' @return A \code{dataframe} maximum genetic distance between samples
#'
#' @rdname dmax


dmax <- function(distMat, popList, pairs, ploidy) {
  if(length(distMat)){
    dmax <- lapply(pairs, function(f){
      #population distance matrix for pairwise pop
      popD <- distMat[as.vector(outer(popList[[f[1]]]$Sample, 1:ploidy, paste, sep = "/")),
                      as.vector(outer(popList[[f[2]]]$Sample, 1:ploidy, paste, sep = "/"))]
      #make a tibble with the average number of pairwise differences
      dmax <- tibble(max(popD))
      #name col
      colnames(dmax) <-  paste0(f[1], ":" , f[2], c("_dmax"))
      dmax
    })

    bind_cols(dmax)
  }
}
