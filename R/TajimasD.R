TajimaD <- function(genoMat, popList, countN = FALSE){

    lapply(popList, function(x){

      smMat <- genoMat[,as.vector(outer(x$Sample, 1:ploidy, paste, sep = "/"))]


      uniq <- apply(smMat , 1, function(x)length(unique(x)))


      # Does the user want missing data to count as differences
      if(!countN){
      Ns <- rowSums(smMat == "N") > 0
      ss <- sum((uniq - Ns) > 1)
      } else{
        ss <- sum(uniq > 1)
      }

      geno <- union(smMat, smMat)
      diff(smMat)

      geno <- geno[geno != "N"]


    dat <- Reduce("+", dat)
    nonN <- rowSums(smMat != "N")

   nonN - dat

    S <- sum(( ncol(genoMat) - dat ) > 0)
  })



}
