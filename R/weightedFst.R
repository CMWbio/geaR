#' Fst from hamming distance
#'
#' @description Calculates sample size weighted or unweighted Fst-like \eqn{\gamma}st from Nei (1982)
#'
#' @details Authours: Chris Ward Uses the mean nucleotide diversity within each subpopulaiton to calculate subpopulaiton pi and divides by the total nucleotide diversity.
#'
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param pairs list of populaiton pairs generated in \code{getDiversityStats}
#' @param ploidy \code{numeric} number of chromosomes
#'
#' @param weighted \code{logical} \cr
#' Weight by sample subpopulaiton size.
#' If \code{TRUE} mean nuceotide diversity within the subpopulaiton will be weigthed by sample size. \cr
#' With weighting 6 samples in p1 and 4 in p2, will result in 0.6*pi_p1 + 0.4*pi_p2 \cr
#' Alternatively, without weighting each subpopulaiton is assumed to contribute equally to the mean, (pi_p1 + pi_p2)/2
#'
#'
#'
#' @return A \code{dataframe} Fst between two populations
#'
#' @rdname Nei82Fst


Nei82Fst <- function(distMat, popList, pairs, ploidy, weighted = TRUE) {
  if(length(distMat)){
    weightedFst <- lapply(pairs, function(f){

      p1 <- popList[[f[1]]]$Sample

      p2 <- popList[[f[2]]]$Sample

      p1alleles <- as.vector(outer(p1, 1:ploidy, paste, sep = "/"))

      p2alleles<- as.vector(outer(p2, 1:ploidy, paste, sep = "/"))

      allP <- c(p1alleles, p2alleles)


      if(weighted) {
        w <- length(p1)/(length(p1)+length(p2))
        piB1 <- w * mean(distMat[p1alleles, p1alleles], na.rm = TRUE)
        piB2 <- (1 - w) * mean(distMat[p2alleles, p2alleles], na.rm = TRUE)

        piB <- piB1 + piB2

      } else {

        piB1 <- mean(distMat[p1alleles, p1alleles], na.rm = TRUE)
        piB2 <- mean(distMat[p2alleles, p2alleles], na.rm = TRUE)

        piB <- piB1 + piB2 / 2
      }

      piT <- mean(distMat[allP, allP], na.rm = TRUE)


      ## as pi
      fst <- tibble(1 - piB / piT)

      #name col
      colnames(fst) <-  paste0(f[1], "_vs_" , f[2], c("_Fst"))
      fst
    })

    bind_cols(weightedFst)
  }
}
