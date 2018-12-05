#' Calculates weighted Fst
#'
#' @description Fst from hamming distance
#'
#' @details Authours: Chris Ward
#' calculates Fst between all populations
#'
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param pairs list of populaiton pairs generated in \code{getDiversityStats}
#' @param ploidy \code{numeric} number of chromosomes
#'
#'
#' @return A \code{dataframe} Fst between two populations
#'
#' @rdname weightedFst


weightedFst <- function(distMat, popList, pairs, ploidy) {
  if(length(distMat)){
    weightedFst <- lapply(pairs, function(f){

      p1 <- popList[[f[1]]]$Sample

      p2 <- popList[[f[2]]]$Sample

      p1alleles <- as.vector(outer(p1, 1:ploidy, paste, sep = "/"))

      p2alleles<- as.vector(outer(p2, 1:ploidy, paste, sep = "/"))

      allP <- c(p1alleles, p2alleles)

      w <- length(p1)/(length(p1)+length(p2))

      piB <- w*mean(distMat[p1alleles, p1alleles]) +
        (1-w)*mean(distMat[p2alleles, p2alleles])


      piT <- mean(distMat[allP, allP])

      fst <- data_frame(1 - piB / piT)

      #name col
      colnames(fst) <-  paste0(f[1], "v" , f[2], c("_Fst"))
      fst
    })

    bind_cols(weightedFst)
  }
}
