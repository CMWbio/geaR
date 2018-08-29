neisPi <- function(dist, popList, ploidy = 2){

  pi <- lapply(popList, function(x){
    #make population matrix
    popD <- distMat[as.vector(outer(x$Sample, 1:ploidy, paste, sep = "/")),
                    as.vector(outer(x$Sample, 1:ploidy, paste, sep = "/"))]
    n <- ncol(popD)
    # get nucleotide diversity
    pi <- sum(popD)/(n*(n-1)/2)
    #determine pi for pop
    piDF <- data_frame(pi)
    # set colnames
    colnames(piDF) <- paste(x$Population[1],"pi", sep = "_")
    piDF
  })

  bind_cols(pi) #put all onto one row

}
