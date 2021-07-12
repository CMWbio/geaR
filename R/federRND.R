#' Calculates Relative Node Depth
#'
#' @description Calculates Relative Node Depth from Feder (2005)
#'
#' @details Authours: Chris Ward
#'
#'
#' @param distMat hamming distance matrix calculated using the hidden function \code{genoDist}
#' @param popList List of populations made from \code{pops} dataframe provided by users
#' @param pairs list of populaiton pairs generated in \code{getDiversityStats}
#' @param ploidy \code{numeric} number of chromosomes
#' @param outgroup 2
#'
#'
#' @return A \code{dataframe} containing relative node depth
#'
#' @rdname RND

RND <- function(distMat, popList, pairs, ploidy, outgroup, type = "feder") {
  if(length(distMat)){
    RND <- lapply(pairs, function(f){

      if(!outgroup %in% f) {


        p1 <- f[1]
        p2 <- f[2]

        #calculate mean distances to outgroup
        dxo <- distMat[as.vector(outer(popList[[p1]]$Sample, 1:ploidy, paste, sep = "/")),
                       as.vector(outer(popList[[outgroup]]$Sample, 1:ploidy, paste, sep = "/"))]

        dyo <- distMat[as.vector(outer(popList[[p2]]$Sample, 1:ploidy, paste, sep = "/")),
                       as.vector(outer(popList[[outgroup]]$Sample, 1:ploidy, paste, sep = "/"))]

        dout <- (mean(dxo) + mean(dyo)) / 2

        if(type == "feder"){

          #calculate dxy
          dxy <- mean(distMat[as.vector(outer(popList[[p1]]$Sample, 1:ploidy, paste, sep = "/")),
                              as.vector(outer(popList[[p2]]$Sample, 1:ploidy, paste, sep = "/"))])


          RND <- tibble::tibble(dxy / dout)

          ## in formula RND(I1,I2;O)
          colnames(RND) <-  paste0("RND(",  p1, ",", p2, ";", outgroup, ")")
          RND

        }

        if(type == "min"){

          #calculate dxy
          dmin <- mean(distMat[as.vector(outer(popList[[p1]]$Sample, 1:ploidy, paste, sep = "/")),
                                as.vector(outer(popList[[p2]]$Sample, 1:ploidy, paste, sep = "/"))])

          RND <- tibble::tibble(dmin / dout)

          ## in formula RND(I1,I2;O)
          colnames(RND) <-  paste0("RNDmin(",  p1, ",", p2, ";", outgroup, ")")
          RND

         }

        RND
      }
    })

   dplyr::bind_cols(RND)
  }

}
