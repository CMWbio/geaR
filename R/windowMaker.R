#' Sliding and tiled window ranges
#'
#' @description makes sliding and tiled windows from data frame
#'
#' @details Authours: Chris Ward & Alastair Ludington
#' generates tiled (no step size) and sliding (step size) ranges for analysis with other functions.
#'
#'
#' @param contigMD data frame with contig metadata. Contains two columns: \code{ID} with scaffold names, \code{length} with scaffold length
#' @param WindowSize \code{numeric} Window size in base pairs
#' @param stepSize \code{numeric} Window step size in base pairs. Default is 0, specifies a tiled window.
#'
#'
#' @return A \code{list} of \code{GRanges}
#'
#'
#' @import SeqArray
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges slidingWindows
#'
#'
#' @examples
#'
#' # Set up tiled windows
#' windowSize <- 100000
#' stepSize <- 0
#'
#' # make contig metadata across all contigs in VCF header
#' #contigMD <- seqVCF_Header(vcf_path)$contig
#'
#' # subset by desired contigs
#' #contigMD <- contigMD[contigMD$ID %in% c("chr1", "chr2"),]
#'
#'
#' @export
#' @rdname windowMaker

windowMaker <- function(contigMD, windowSize, stepSize = 0, nCores = 1){

  # checks
  if(!is.data.frame(contigMD)) stop("contigMD must be a data.frame or data.frame like object")
  if(any(!colnames(contigMD) %in% c("ID", "length"))) stop("contigMD must be a data.frame or data.frame like object")
  if(!windowSize >= stepSize) stop("windowSize must be greater than or equal to stepSize")
  if(!is.numeric(windowSize) | !is.numeric(stepSize)) stop("windowSize and stepSize must be numerics")
  stopifnot(windowSize > 0 | stepSize >= 0 )

  # set step size for tiled
  if(stepSize == 0) stepSize <- windowSize

  # remove all contigs with length less than windowSize
  contigMD <- contigMD[contigMD[["length"]] >= windowSize,]

  # make contig range for scaffolds in contigMD and convert to GRanges object
  contigRange <- makeGRangesFromDataFrame(data_frame(chr = contigMD[["ID"]], start = 1, strand = ".", end = contigMD[["length"]]))


  ## List of windows
  windowList <- mclapply(seq(nrow(contigMD)), mc.cores = nCores, function(x){

    ## Subsetting GRange object by chromosome
    con <- contigRange[x]

    ## Splitting into windows
    G_range <- GenomicRanges::slidingWindows(x = con, width = windowSize, step = stepSize)
    G_range <- unlist(G_range) ## Need to unlist the slidingWindow object

    ## Turning window-GRanges into GRangesList

  })

  ## this is slow
  windowList <- do.call("c",windowList)
  windowList <- as(windowList,"GRangesList")
  return(windowList)

}


