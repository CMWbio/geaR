TajimaD <- function(genoMat, popList, countN){

  lapply(popList, function(x){

    geno <- union(genoMat, genoMat)

    if(!countN) geno <- geno[geno != "N"]

    dat <- lapply(geno, function(f){

      mat <- genoMat == f
      rowsum()
    })

    dat <- Reduce("+", dat)

    distMat <- ( nrow(genoMat) - dat ) / nrow(genoMat)

  })



}
