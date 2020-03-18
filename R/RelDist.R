#' Calculates Relative genetic distance (dxy)
#'
#' @description Calculates Relative distance, an extension of ttd from Ward and Baxter 2018
#'
#' @details 
#'
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param ploidy \code{numeric} number of chromosomes
#' @param outgroup \code{character} the furthest outgroup
#'
#'
#' @return A \code{dataframe} containing Relative distance
#'
#' @rdname relDxy

relDxy <- function(distMat, popList, ploidy, outgroup){
    
    pops <- bind_rows(popList)
    uTest <- unique(pops$Population)
    noOut <- uTest[uTest != outgroup]
    trips <- combinat::permn(noOut)
    
    allD <- lapply(trips, function(x){
        
        di1o1 <- mean(distMat[as.vector(outer(popList[[x[1]]]$Sample, 1:ploidy, paste, sep = "/")),
                       as.vector(outer(popList[[x[3]]]$Sample, 1:ploidy, paste, sep = "/"))])
        
        di2o1 <- mean(distMat[as.vector(outer(popList[[x[2]]]$Sample, 1:ploidy, paste, sep = "/")),
                       as.vector(outer(popList[[x[3]]]$Sample, 1:ploidy, paste, sep = "/"))])
        
        di1o2 <- mean(distMat[as.vector(outer(popList[[x[1]]]$Sample, 1:ploidy, paste, sep = "/")),
                         as.vector(outer(popList[[outgroup]]$Sample, 1:ploidy, paste, sep = "/"))])
        
        di2o2 <- mean(distMat[as.vector(outer(popList[[x[2]]]$Sample, 1:ploidy, paste, sep = "/")),
                         as.vector(outer(popList[[outgroup]]$Sample, 1:ploidy, paste, sep = "/"))])
        
        
        rD <- (di2o1/(di2o1+di1o1))/(di2o2/(di2o2+di1o2))
        
        
        df <- tibble(rD)
        colnames(df) <- paste0("rDxy(", x[1], ",", x[2], ",", x[3], ";", outgroup, ")")
        df
    
    })
    
    bind_cols(allD)
}