#' Calculates minimum genetic distance
#'
#' @description Calculates minimum population genetic distance from hamming distance
#'
#' @details Authours: Chris Ward
#' Calculates minimum population genetic distance from hamming distance
#'
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param pairs list of populaiton pairs generated in \code{getDiversityStats}
#' @param ploidy \code{numeric} number of chromosomes
#'
#'
#' @return A \code{dataframe} mean minimum genetic distance between samples
#'
#' @rdname dmin
#' @export


dmin <- function(distMat, popList, pairs, ploidy) {
  if(length(distMat)){
    dmin <- lapply(pairs, function(f){
      #population distance matrix for pairwise pop
      popD <- distMat[as.vector(outer(popList[[f[1]]]$Sample, 1:ploidy, paste, sep = "/")),
                      as.vector(outer(popList[[f[2]]]$Sample, 1:ploidy, paste, sep = "/"))]
      #make a tibble with the average number of pairwise differences
      dmin <- data_frame(min(popD))
      #name col
      colnames(dmin) <-  paste0(f[1], "v" , f[2], c("_dmin"))
      dmin
    })

    bind_cols(dmin)
  }
}
