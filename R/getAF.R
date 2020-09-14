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

getAF <- function(GDS, locus, minSites, pops, refAllele, counts = FALSE){
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

    seqname <- seqGetData(GDS, var.name = "chromosome")
    ## Filtering empty arrays/arrays with too few variants
    if(length(pos) < minSites) {

      print("Incompatible window - not enough information")
      return(NULL)
    }
    else{



      ref <- seqGetData(GDS, var.name = "$ref")
      alt <- seqGetData(GDS, var.name = "$alt")

      if(counts) AF <- seqAlleleCount(GDS, ref.allele =  refAllele)
      else AF <- seqAlleleFreq(GDS, ref.allele =  refAllele)
      

      df <- tibble(seqname,pos,ref,alt,AF)
      df$alt[df$alt == ""] <- "N"
      df <- df[nchar(df$ref) == nchar(df$alt) & nchar(df$alt) == 1,]
      colnames(df) <- c("seqname", "pos", "ref", "alt", paste0(n, "_AF"))

      df


    }

  })

  data <- Filter(Negate(is.null), data)
  if(length(data)) data <- reduce(data, left_join, by = c("seqname","pos", "ref", "alt"))



}
