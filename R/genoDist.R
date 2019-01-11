#' Calculates the hamming distance between alleles
#'
#' @description Used in the calculation of diversity statistics
#'
#' @details Authours: Chris Ward
#' Calculates the hamming distance using matrix multiplication
#'
#'
#' @param genoMat A \code{matrix} \cr
#' Allele genotypes for each individual
#' @param pairwiseDeletion  \code{logical}
#' If \code{TRUE} missing data will be removed from distance calculations in a paiwise manner.
#'
#' @return A \code{matrix} of hamming distance between individuals
#'
#' @rdname genoDist
#' @export

genoDist <- function(genoMat, pairwiseDeletion){

  geno <- union(genoMat, genoMat)
  dat <- lapply(geno, function(f){
    mat <- genoMat == f
    t(mat) %*% mat
  })

  names(dat) <- geno


  sim <- Reduce("+", dat)

  S <- nrow(genoMat)
  dif <- S - sim


  if(pairwiseDeletion & any(genoMat == "N")){

  notNmat <- genoMat != "N"
  ## get number of non N sites in each pairwise
  nonNsites <- t(notNmat) %*% notNmat
  #total number of Ns between each individual
  Ntotal <-  S - nonNsites
  # N same between pairwise indv
  Nsim <- dat[["N"]]
  # N different between pairwise indv
  Ndif <- Ntotal - Nsim
  # non N differences
  nonNdif <- dif - Ndif
  # differences without Ns
  distMat <- nonNdif / nonNsites


  }else{
  distMat <- dif / nrow(genoMat)

    }

  return(distMat)



}
