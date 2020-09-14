#' Gmin
#'
#' @description Calculates Gmin from Geneva (2015)
#'
#' @details Authours: Chris Ward \cr
#' calculates Gmin from ...
#'
#'
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param pairs list of populaiton pairs generated in \code{getDiversityStats}
#' @param ploidy \code{numeric} number of chromosomes
#'
#'
#' @return A \code{dataframe} containing Gmin
#'
#' @rdname Gmin
#' @export


Gmin <- function(distMat, popList, pairs, ploidy) {
  if(length(distMat)){
    Gmin <- lapply(pairs, function(f){

      p1 <- f[1]
      p2 <- f[2]

        #calculate d
        d <- distMat[as.vector(outer(popList[[p1]]$Sample, 1:ploidy, paste, sep = "/")),
                            as.vector(outer(popList[[p2]]$Sample, 1:ploidy, paste, sep = "/"))]
        #calculate dxy and dmin
        dxy <- mean(d)

        dmin <- min(d)

        Gmin <- tibble(dmin / dxy)

        ## in formula RND(I1,I2;O)
        colnames(Gmin) <-  paste0(p1, "v", p2, "_Gmin")
        Gmin
    })

    dplyr::bind_cols(Gmin)
  }

}
