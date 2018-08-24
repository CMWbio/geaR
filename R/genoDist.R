genoDist <- function(genoMat){

  geno <- union(genoMat, genoMat)

  dat <- lapply(geno, function(f){

    mat <- genoMat == f
    t(mat) %*% mat
  })

  dat <- Reduce("+", dat)

  distMat <- ( nrow(genoMat) - dat ) / nrow(genoMat)

}


