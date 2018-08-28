neisDxy <- function(distMat, popList, pairs, ploidy) {
  if(length(distMat)){
    dxy <- lapply(pairs, function(f){
          #population distance matrix for pairwise pop
          popD <- distMat[as.vector(outer(popList[[f[1]]]$Sample, 1:ploidy, paste, sep = "/")),
                          as.vector(outer(popList[[f[2]]]$Sample, 1:ploidy, paste, sep = "/"))]
          #make a tibble with the average number of pairwise differences
          dxy <- data_frame(mean(popD), sd(popD))
          #name col
          colnames(dxy) <-  paste0(f[1], "v" , f[2], c("_dxy", "_SDdxy"))
          dxy
    })

    bind_cols(dxy)
  }
}
