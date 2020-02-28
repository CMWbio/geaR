#' Sliding or Tiled snpwindow maker.
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
#' @import tibble
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges slidingWindows
#' @importFrom parallel mclapply
#'
#' @name snpWindowMaker
#' @rdname snpWindowMaker-methods
#' @export

setGeneric("snpWindowMaker", function(x, windowSize, stepSize, nCores = 1, GDS = NULL, snpWindow = FALSE, loci = NULL, ...){standardGeneric("snpWindowMaker")})

#' @aliases snpWindowMaker,data_frame
#' @rdname snpWindowMaker-methods
#' @export
setMethod("snpWindowMaker", signature = "data.frame",
          function(x, windowSize, stepSize, nCores, GDS = NULL, snpWindow = FALSE, loci = NULL){

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
            windowList <- mclapply(seq(length(contigRange)), mc.cores = nCores, function(y){

              ## Subsetting GRange object by chromosome
              con <- contigRange[y]

                seqSetFilter(object = GDS, variant.sel = con)
                position <- seqGetData(gdsfile = GDS, var.name = "position")
                gr <- GRanges(seqnames = con@seqnames, IRanges(start = position, end = position))

                if(!is.null(loci)){

                  gr <- join_overlap_inner(gr, loci@unlistData)

                }


                position <- gr@ranges@start
                if(length(position) > windowSize) {
                  position <- zoo::rollapply(position, width = windowSize, by = stepSize, function(x){
                    x
                  })

                  posList <- split(position, 1:nrow(position))


                  grL <- lapply(posList, function(posX){

                    gr <- gr[gr@ranges@start %in% posX]

                    if(stepSize == windowSize) gr$lociType<- "tiledSnpWindow"
                    else gr$lociType<- "slidingSnpWindow"

                    gr
                  })

                  grL


                }




            })

            windowList <- do.call("c", windowList)

            GRangesList(windowList)

          }
)
