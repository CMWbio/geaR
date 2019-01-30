#' Combines overlaps in features and windows
#'
#' @description Combines multiple GRange objects with overlapping ranges.
#'
#' @details Authour: Chris Ward
#' Passing multiple \code{GRange} objects will result in ordered overlaping of ranges. e.g. if a genome wide 100kb tiled window \code{GRangeList}
#'  generated using \code{windowMaker()} is passed along with a features \code{GRangesList} generated using \code{getFeatures(feature = "gene:cds")},
#'  This will be combined into 100kb windows containing only regions of protein coding sequence.
#'
#'
#' @param window  \code{GRangesList}. Windows to overlap features to
#' @param feature \code{GRangesList} containing features generated using either \code{getFeatures} or \code{codonRanges}
#' @param nCores \code{numeric} number of cores to run in parallel
#' @param overlap \code{character} \cr
#' Negative or positive overlap
#' "-" will subtract the feature ranges from the windows \cr
#' "+" will bin features into window ranges
#'
#'
#' @return A \code{data_frame} of selected Diversity statistics
#'
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#'
#' @export
#' @rdname combineRanges


combineRanges <- function(window, feature, nCores, overlap = "+"){

    data <- mclapply(window, mc.cores = nCores, function(x){

    overlapSeqName <- feature@unlistData[as.character(feature@unlistData@seqnames) == as.character(x@seqnames)]

    grangeInRange <- overlapSeqName[overlapSeqName@ranges@start >= x@ranges@start & end(overlapSeqName) <= end(x),]

    if(!length(grangeInRange)) grangeInRange <- NA

    if(overlap == "-" & any(!is.na(grangeInRange))) {

      grangeInRange  <- psetdiff(x, GRangesList(grangeInRange))
      grangeInRange <- unlist(grangeInRange)


    }

    grangeInRange

    })



  data <- data[!is.na(data)]

  GRangesList(data)

}
