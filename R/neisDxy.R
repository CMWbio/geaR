neisDxy <- function(distMat, popList, ploidy) {
  if(length(distMat)){
  dxy <- lapply(popList, function(x){
    dInner <-lapply(popList, function(y){
      if(!setequal(x$Sample, y$Sample)){
        #population distance matrix for pairwise pop
        popD <- distMat[as.vector(outer(as.character(x$Sample), 1:ploidy, paste, sep = "/")),
                        as.vector(outer(as.character(y$Sample), 1:ploidy, paste, sep = "/"))]
        #make a tibble with the average number of pairwise differences
        dxy <- data_frame(mean(popD), sd(popD))
        #name col
        colnames(dxy) <-  paste0(x$Population[1], "v/" , y$Population[1], c("_dxy", "_SDdxy"))
        dxy
      }
    })
    bind_cols(dInner)
  })

  bind_cols(dxy)
  }
}
