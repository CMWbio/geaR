#' Calculates population based allele frequency
#'
#' @description Calculates population based allele frequency
#'
#' @details Authours: Chris Ward
#' Calculates population based allele frequency
#'
#'
#' @param GDS allele frequency
#' @param locus Granges locus to act on
#' @param minSites minimum sites within a window to act on
#' @param pop populations dataframe
#'
#' @return A \code{tibble} of population allele frequency for biallelic SNPs
#'
#' @importFrom purrr reduce
#'
#' @rdname getAF
#' @export

getAF <- function(GDS, locus, minSites, pops){
  stopifnot(minSites < 1 & minSites > 0)
  minSites <- sum(minSites * width(locus))

  popL <- split(pops, pops$Population)

  data <- map(seq(length(popL)), function(x){

    n <- names(popL)[x]
    x <- popL[[x]]

    samples <- x$Sample
    seqSetFilter(object = GDS, sample.id = samples)
    seqSetFilter(object = GDS, variant.sel = locus)

    pos <- seqGetData(GDS, var.name = "position")

    ## Filtering empty arrays/arrays with too few variants
    if(length(pos) < minSites) {

      print("Incompatible window - not enough information")
      return(NULL)
    }
    else{



      ref <- seqGetData(GDS, var.name = "$ref")
      alt <- seqGetData(GDS, var.name = "$alt")

      AF <- seqAlleleFreq(GDS)

      df <- tibble::tibble(pos,ref,alt,AF)
      df <- df[nchar(df$ref) == nchar(df$alt) & nchar(df$alt) == 1,]
      colnames(df) <- c("pos", "ref", "alt", paste0(n, "_AF"))

      df


    }

  })

  data <- reduce(data, left_join, by = c("pos", "ref", "alt"))



}
