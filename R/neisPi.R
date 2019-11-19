#' Calculates Nei's nucleotide diversity (pi)
#'
#' @description Calculates Nei (1987) pi from hamming distance
#'
#' @details Authours: Chris Ward
#' calculates within population nucleotide diversity
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param ploidy \code{numeric} number of chromosomes
#'
#'
#' @return A \code{dataframe} mean nucleotide diversity within each populaiton
#'
#' @rdname neisPi
#' @export


neisPi <- function(distMat, popList, ploidy = 2){

  pi <- lapply(popList, function(x){
    #make population matrix
    popD <- distMat[as.vector(outer(x$Sample, 1:ploidy, paste, sep = "/")),
                    as.vector(outer(x$Sample, 1:ploidy, paste, sep = "/"))]
    n <- ncol(popD)
    # get nucleotide diversity

    # if(sampleCorrect) pi <- sum(popD)/(n*(n-1))
    # else
    pi <- mean(popD)
    #determine pi for pop
    piDF <- data_frame(pi)
    # set colnames
    colnames(piDF) <- paste(x$Population[1],"pi", sep = "_")
    piDF
  })

  bind_cols(pi) #put all onto one row

}
