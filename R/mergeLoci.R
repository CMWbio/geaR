#' Combines overlaps in ys and xs
#'
#' @description Combines multiple GRangeList objects with overlapping ranges.
#'
#' @details Authour: Chris Ward
#' Passing multiple \code{GRange} objects will result in ordered overlaping of ranges. e.g. if a genome wide 100kb tiled x \code{GRangeList}
#'  generated using \code{xMaker()} is passed along with a ys \code{GRangesList} generated using \code{getwindows(y = "gene:cds")},
#'  This will be combined into 100kb xs containing only regions of protein coding sequence.
#'
#' @param x \code{GRangesList}. xs to overlap ys to.
#' @param y \code{GRangesList} or \code{character} containing ys. Provide wither a GRangesList or path to SQL database.
#' @param nCores \code{numeric} number of cores to run in parallel
#' @param overlap \code{character} \cr
#' Negative or positive overlap
#' "-" will subtract the y ranges from the xs \cr
#' "+" will bin ys into x ranges
#' @param table \code{character} one of \code{c("4Degen", "0Degen", "2Degen")} 
#' 
#' @return A \code{GRangesList} containing merged ranges 
#'
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom plyranges join_overlap_inner
#' @importFrom plyranges setdiff_ranges
#'
#' @export
#' @rdname mergeLoci


setGeneric("mergeLoci", function(x, y, nCores, overlap = "+", table = NULL){
  standardGeneric("mergeLoci")
})


setMethod("mergeLoci", signature(c(y = "character")),
          function(x, y, nCores, overlap, table){
            
            conn <- RSQLite::dbConnect(RSQLite::SQLite(), y)
            
            if(!is.null(table)){
              
              df <-  RSQLite::dbGetQuery(conn, paste0("SELECT * FROM '", table, "'"))
              
            }
            
            else {
              
              df <-  RSQLite::dbGetQuery(conn, paste0("SELECT * FROM genes"))[[1]]
              
              df <- purrr::map(df, function(x){
                df <-  RSQLite::dbGetQuery(conn, paste0("SELECT seqnames,strand,start,end FROM '", x, 
                                                        "' WHERE codonPosition IN ('", paste(positions, collapse = "', '"),"')"))
              }) 
              
              df <- bind_rows(df)
            }
            
            data <- mclapply(seq(x), mc.cores = nCores, function(locus){
              
              locus <- x[[locus]]
              
              overlapSeqName <- df[as.character(df$seqnames) == as.character(locus@seqnames),]
              
              grangeInRange <- overlapSeqName[overlapSeqName$start >= locus@ranges@start[1] & overlapSeqName$end <= end(locus)[length(locus)],]
              grangeInRange <- makeGRangesFromDataFrame(grangeInRange)
              if(!length(grangeInRange)) grangeInRange <- NA
              
              if(any(!is.na(grangeInRange))) {
                
                if(overlap == "+") join_overlap_inner(grangeInRange, locus )
                else setdiff_ranges(locus, grangeInRange)
                
              }
              
            })
            
            data <- data[!is.na(data)]
            data <- Filter(Negate(is.null), data)
            
            GRangesList(data)
            
          })

setMethod("mergeLoci", signature(c(y = "GRangesList")),
          function(x, y, nCores, overlap, table){
            
            data <- mclapply(seq(x), mc.cores = nCores, function(locus){
              locus <- x[[locus]]
              
              y <- y@unlistData #[y@unlistData$codonPosition %in% positions,]
              y@elementMetadata <- y@elementMetadata[colnames(y@elementMetadata) == "gene"]
              overlapSeqName <- y[as.character(y@seqnames) == as.character(locus@seqnames)[1]]
              
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
            
            
            
          })
