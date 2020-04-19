#' #' Calculates geentic distance between populations using allele frequency at biallelic SNPS
#' #'
#' #' @description Used in the calculation of diversity statistics FAST
#' #'
#' #' @details Authours: Chris Ward
#' #' Calculates geentic distance between populations using allele frequency at biallelic SNPS
#' #'
#' #'
#' #' @param AF allele frequency
#' #' @param pop populaitons dataframe
#' #'
#' #' @return A \code{matrix} of genetic distance between populaitons
#' #'
#' #' @useDynLib geaR
#' #'
#' #'
#' #' @rdname dxyAF
#' #' @export
#'
dxyAF <- function(AF){

    mat <- as.matrix(AF[5:ncol(AF)])
    

    t <- lapply(1:nrow(mat), function(x){
        dist(mat[x,])
        
    })
    
    Reduce('+', t) / nrow(AF)

}
