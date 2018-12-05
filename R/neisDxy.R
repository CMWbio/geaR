#' Calculates Nei's Dxy
#'
#' @description Calculates Nei (1987) Dxy from hamming distance
#'
#' @details Authours: Chris Ward
#' calculates mean absolute genetic distance between two populations
#'
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param pairs list of populaiton pairs generated in \code{getDiversityStats}
#' @param ploidy \code{numeric} number of chromosomes
#'
#'
#' @return A \code{dataframe} mean absolute genetic distance between samples
#'
#' @rdname neisDxy


neisDxy <- function(distMat, popList, pairs, ploidy) {
  if(length(distMat)){
    dxy <- lapply(pairs, function(f){
          #population distance matrix for pairwise pop
          popD <- distMat[as.vector(outer(popList[[f[1]]]$Sample, 1:ploidy, paste, sep = "/")),
                          as.vector(outer(popList[[f[2]]]$Sample, 1:ploidy, paste, sep = "/"))]
          #make a tibble with the average number of pairwise differences
          dxy <- data_frame(mean(popD), sd(popD))
          #name col
          colnames(dxy) <-  paste0(f[1], "v" , f[2], c("_dxy", "_SDdxy"))
          dxy
    })

    bind_cols(dxy)
  }

}
