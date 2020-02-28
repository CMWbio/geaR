#' Combines overlaps in ys and xs
#'
#' @description Combines multiple GRange objects with overlapping ranges.
#'
#' @details Authour: Chris Ward
#' Passing multiple \code{GRange} objects will result in ordered overlaping of ranges. e.g. if a genome wide 100kb tiled x \code{GRangeList}
#'  generated using \code{xMaker()} is passed along with a ys \code{GRangesList} generated using \code{getwindows(y = "gene:cds")},
#'  This will be combined into 100kb xs containing only regions of protein coding sequence.
#'
#' @param x \code{GRangesList}. xs to overlap ys to
#' @param y \code{GRangesList} containing ys generated using either \code{getfeatures} or \code{codonRanges}
#' @param nCores \code{numeric} number of cores to run in parallel
#' @param overlap \code{character} \cr
#' Negative or positive overlap
#' "-" will subtract the y ranges from the xs \cr
#' "+" will bin ys into x ranges
#'
#' @return A \code{data_frame} of selected Diversity statistics
#'
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#'
#' @export
#' @rdname combineRanges


combineRanges <- function(x, y, nCores, overlap = "+"){

    data <- mclapply(x, mc.cores = nCores, function(locus){

    overlapSeqName <- y@unlistData[as.character(y@unlistData@seqnames) == as.character(locus@seqnames)[1]]

    grangeInRange <- overlapSeqName[overlapSeqName@ranges@start >= locus@ranges@start[1] & end(overlapSeqName) <= end(locus)[length(locus)],]

    if(!length(grangeInRange)) grangeInRange <- NA

    if(any(!is.na(grangeInRange))) {

     if(overlap == "+") join_overlap_inner(grangeInRange, x )
     else setdiff_ranges(x, grangeInRange)

    }

    })



  data <- data[!is.na(data)]
  data <- Filter(Negate(is.null), data)

  GRangesList(data)

}
