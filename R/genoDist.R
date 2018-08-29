genoDist <- function(genoMat){

  geno <- union(genoMat, genoMat)

  dat <- lapply(geno, function(f){

    mat <- genoMat == f
    t(mat) %*% mat
  })



  sim <- Reduce("+", dat)

  dif <- nrow(genoMat) - sim

  distMat <- dif / nrow(genoMat)



  ## started to think about how to remove Ns, will have to itterate through and pairwise delete
  ## I am not smart enought to figure out a matrix solution
  # if(!countN){
  # ## get number of N between each
  # Ns <- t(genoMat != "N" ) %*% (genoMat != "N")
  # Ntotal <-  nrow(genoMat) - Ns
  #
  # # correct using proportions
  # Nprop <- Ntotal / nrow(genoMat)
  #
  # distMat - Npr
  #
  #
  # # # get N corrected dif and effective number of sites
  # # dif <- dif - Ntotal
  # # #dif[dif < 0] <- 0
  # #
  # # sitesEff <- nrow(genoMat) - Ntotal
  # #
  # # distMat <- dif / sitesEff
  #
  # }else{
  #
  #
  #
  #   }
  #
  return(distMat)



}


