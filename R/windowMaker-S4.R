#' Sliding or Tiled window maker.
#'
#' @description Makes sliding and tiled windows from data frame
#'
#' @details Authours: Chris Ward & Alastair Ludington \cr
#' Generates tiled (no step size) and sliding (step size) ranges for analysis with other functions.
#'
#'
#' @param x Can be a \code{data.frame} object containing ID and Length columns for each scaffold, \cr
#'  \code{Granges} or \code{GRangesList} object
#' @param WindowSize \code{numeric} \cr Window size in base pairs
#' @param stepSize \code{numeric} \cr Window step size in base pairs. Default is 0, specifies a tiled window.
#'
#'
#' @return A \code{GRangesList} object where each element is a genomic window.
#'
#'
#' @import SeqArray
#' @importFrom tibble data_frame
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges slidingWindows
#' @importFrom parallel mclapply
#'
#'
#' @examples
#'
#' #### tiled windows ####
#' # make tiled windows of width 10,000 bp using 4 cores
#' windowSize <- 10000
#' stepSize <- 0
#' nCores <- 4
#'
#' #### sliding windows ####
#' # make sliding windows of width 10,000 bp with 1000 bp step size,  using 4 cores
#' windowSize <- 10000
#' stepSize <- 1000
#' nCores <- 4
#'
#'
#' #### Using data.frame ####
#'
#' # get contig metadata from fasta index
#' df <- readr::read_tsv("referene.fna.fai", col_names = FALSE)
#' x <- tibble::data_frame(ID = df[[1]], length = df[[2]])
#'
#' # subset by desired contigs
#'
#' x <- x[x$ID %in% c("chr1", "chr2"),]
#'
#' # construct windows
#'
#' windows <- windowMaker(x, windowSize, stepSize, nCores)
#'
#' #### Using GRanges ####
#'
#' # make GRange object of width 1 Mbp on chr1
#'
#' x <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 1000000))
#' windows <- windowMaker(x, windowSize, stepSize)
#'
#' #### Using GRangesList ####
#'
#' # make GRange object of width 1 Mbp on chr1 and chr2
#' x <-GRangesList(GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 1000000)), GRanges(seqnames = "chr2", ranges = IRanges(start = 1, end = 1000000)))
#'
#' windows <- windowMaker(x, windowSize, stepSize, nCores)
#'
#' @name windowMaker
#' @rdname windowMaker-methods
#' @export
#'

setGeneric("windowMaker",function(x, windowSize, stepSize, ...){standardGeneric("windowMaker")})

#' @aliases windowMaker,data_frame
#' @rdname windowMaker-methods
#' @export
setMethod("windowMaker", signature = "data.frame",
          function(x, windowSize, stepSize = 0, nCores = 1){

            # checks
            if(!is.data.frame(x)) stop("x must be a data.frame or data.frame like object")
            if(any(!colnames(x) %in% c("ID", "length"))) stop("x must be a data.frame or data.frame like object")
            if(!windowSize >= stepSize) stop("windowSize must be greater than or equal to stepSize")
            if(!is.numeric(windowSize) | !is.numeric(stepSize)) stop("windowSize and stepSize must be numerics")
            stopifnot(windowSize > 0 | stepSize >= 0 )

            # set step size for tiled
            if(stepSize == 0) stepSize <- windowSize

            # remove all contigs with length less than windowSize
            x <- x[x[["length"]] >= windowSize,]

            # make contig range for scaffolds in x and convert to GRanges object
            contigRange <- makeGRangesFromDataFrame(data_frame(chr = x[["ID"]], start = 1, strand = ".", end = x[["length"]]))


            ## List of windows
            windowList <- mclapply(seq(nrow(x)), mc.cores = nCores, function(y){

              ## Subsetting GRange object by chromosome
              con <- contigRange[y]

              ## Splitting into windows
              G_range <- slidingWindows(x = con, width = windowSize, step = stepSize)
              G_range <- unlist(G_range) ## Need to unlist the slidingWindow object

              elementMetadata(G_range)[["lociType"]] <- "window"

              G_range



            })

            windowList <- do.call("c", windowList)

            windowList <- split(windowList, as.factor(windowList))

            return(windowList)

          }
)

#' @aliases windowMaker,GRanges
#' @rdname windowMaker-methods
#' @export
setMethod("windowMaker", signature = "GRanges",
          function(x, windowSize, stepSize = 0){

            if(!windowSize >= stepSize) stop("windowSize must be greater than or equal to stepSize")
            if(!is.numeric(windowSize) | !is.numeric(stepSize)) stop("windowSize and stepSize must be numerics")
            if(stepSize == 0) stepSize <- windowSize

            G_range <- slidingWindows(x = x, width = windowSize, step = stepSize)
            G_range <- unlist(G_range)

            elementMetadata(G_range)[["lociType"]] <- "window"

            windowList <- split(G_range, as.factor(G_range))

            return(windowList)

            }
)

#' @aliases windowMaker,GRangesList
#' @rdname windowMaker-methods
#' @export
setMethod("windowMaker", signature = "GRangesList",
          function(x, windowSize, stepSize = 0, nCores = 1){

            if(!windowSize >= stepSize) stop("windowSize must be greater than or equal to stepSize")
            if(!is.numeric(windowSize) | !is.numeric(stepSize)) stop("windowSize and stepSize must be numerics")
            if(stepSize == 0) stepSize <- windowSize

            windowList <- mclapply(seq(length(x)), mc.cores = nCores, function(y){

            G_range <- x[[y]]

            G_range <- slidingWindows(x = G_range, width = windowSize, step = stepSize)
            G_range <- unlist(G_range)

            elementMetadata(G_range)[["lociType"]] <- "window"
            G_range


            })

            windowList <- do.call("c", windowList)

            windowList <- split(windowList, as.factor(windowList))

            return(windowList)

          }
)

