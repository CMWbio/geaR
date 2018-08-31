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
#'

genoDist <- function(genoMat, pairwiseDeletion){

  geno <- union(genoMat, genoMat)
  dat <- lapply(geno, function(f){
    mat <- genoMat == f
    t(mat) %*% mat
  })

  sim <- Reduce("+", dat)

  dif <- nrow(genoMat) - sim




  # started to think about how to remove Ns, will have to itterate through and pairwise delete
  # I am not smart enought to figure out a matrix solution
  if(!countN){
  ## get number of non N sites in each pairwise
  nonNsites <- t(genoMat != "N" ) %*% (genoMat != "N")
  #total number of Ns between each individual
  Ntotal <-  nrow(genoMat) - nonNsites
  # N same between pairwise indv
  N <- genoMat == "N"
  Nsim <- t(N) %*% N
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

#' Calculates teh hamming distance between alleles
#'
#' @description Used in the calculation of diversity statistics
#'
#' @details Authours: Chris Ward
#' Calculates the hamming distance using matrix multiplication
#'
#'
#' @param genoMat A \code{matrix} containing alleles for each individual
#'
#' @return A \code{matrix} of hamming distance between individuals
#'
#'

genoDist2 <- function(genoMat){

  geno <- union(genoMat, genoMat)

  nam <- pops$Sample
  pairs <- outer(nam, nam, paste, sep = "///")
  pairs <- lapply(pairs, function(z){
    pair <- unlist(strsplit(z, split = "///"))
    if(pair[1] == pair[2]) pair <- NULL
    pair <- sort(pair)
    return(pair)
  })
  pairs <- pairs[!duplicated(pairs)]
  pairs <- Filter(Negate(is.null), pairs)

  lapply()


  return(distMat)



}

