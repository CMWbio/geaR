#' Calculates teh hamming distance between alleles
#'
#' @description Used in the calculation of diversity statistics
#'
#' @details Authours: Chris Ward
#' Calculates the hamming distance using matrix multiplication
#'
#'
#' @param genoMat A \code{matrix} containing alleles for each individual
#' @param pairwiseDeletion  \code{logical}
#'
#' @return A \code{matrix} of hamming distance between individuals
#'
#' @rdname genoDist

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




  # started to think about how to remove Ns, will have to itterate through and pairwise delete
  # I am not smart enought to figure out a matrix solution
  if(pairwiseDeletion){

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


  # # get N corrected dif and effective number of sites
  # dif <- dif - Ntotal
  # #dif[dif < 0] <- 0
  #
  # sitesEff <- nrow(genoMat) - Ntotal
  #
  # distMat <- dif / sitesEff

  }else{
  distMat <- dif / nrow(genoMat)

    }

  return(distMat)



}
