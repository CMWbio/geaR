#' filter codon position and type
#'
#'
#' @description 
#' Use a codon database to filter codons based on position and degeneracy
#'
#' @param DB \code{GRangeList} or \code{character} vector containing the path to a SQL database
#' @param position \code{character} any combination of \code{c("first", "second", "third")}.
#' @param degeneracy \code{character} this will select positions based on degeneracy can be one of \code{c(0, 2)}. \cr 
#' This will overwrite the \code{position} argument. \cr 
#' For four fold degenerate sites use \code{valicate4fold()}.  
#'
#'
#' @return An object of class \code{GRanges} containing codon ranges for positions specified.
#'
#'
#' @examples
#'
#' @export
#' @rdname filterCodonDB

setGeneric("filterCodonDB", function(DB, position, degeneracy){
    standardGeneric("filterCodonDB")
})

#' @aliases filterCodonDB,GrangesList
#' @export
setMethod("filterCodonDB", signature(DB = "GRangesList"),
          function(DB, position, degeneracy = NULL){
              
              if(is.null(degeneracy)) return(DB@unlistData[DB@unlistData$codonPosition %in% position])
              degeneracy <- as.numeric(degeneracy)
              
              if(degeneracy == 0){
                  position <- c("first", "second")
                  filt <- c("Trp", "Met")
                  
                  DB@unlistData <- DB@unlistData[DB@unlistData$codonPosition %in% position |
                                                     DB@unlistData$residue %in% filt]
                  
                  DB@unlistData <- DB@unlistData[DB@unlistData$residue != "Stp"]
                  return(DB)
              } 
              
              if(degeneracy == 2){
                  
                  position <- "third"
                  filt <- c("Phe", "Leu2", "Tyr", "His", "Glu", "Asn", "Lys", "Asp", "Glu", "Cys", 
                            "Ser", "Arg2")
                  
                  DB@unlistData <- DB@unlistData[DB@unlistData$residue %in% filt]
                  
                  DB@unlistData <- DB@unlistData[DB@unlistData$codonPosition %in% position]
                  return(DB)
              }
              
              
              
          })

#' @aliases filterCodonDB,character
#' @export
setMethod("filterCodonDB", signature(DB = "character"),
          function(DB, position, degeneracy = NULL){
              
              name <- DB
              conn <- dbConnect(SQLite(), DB)
              genes <- dbGetQuery(conn, paste0("SELECT * FROM genes"))[[1]]
              
              DB <- mclapply(genes, mc.cores = nCores, function(x){
                  # x <- genes[[x]]
                  
                  gene <- dbGetQuery(conn, paste0("SELECT * FROM '", x, "'"))
                  
              })
              
              DB <- dplyr::bind_rows(DB)
              
              
              if(is.null(degeneracy)) {
                  DB <- DB[DB$codonPosition %in% position,]
                  dbWriteTable(conn = conn, name = paste(position, collapse = "_"), value = data, overwrite = TRUE)
                  
              }
              
              if(degeneracy == 0){
                  position <- c("first", "second")
                  filt <- c("Trp", "Met")
                  
                  DB <- DB[DB$codonPosition %in% position |
                               DB$residue %in% filt,]
                  
                  DB <- DB[DB$residue != "Stp",]
                  
                  dbWriteTable(conn = conn, name = "0Fold", value = data, overwrite = TRUE)
                  
              } 
              
              if(degeneracy == 2){
                  
                  position <- "third"
                  filt <- c("Phe", "Leu2", "Tyr", "His", "Glu", "Asn", "Lys", "Asp", "Glu", "Cys", 
                            "Ser", "Arg2")
                  
                  DB <- DB[DB$residue %in% filt,]
                  
                  DB <- DB[DB$codonPosition %in% position,]
                  
                  dbWriteTable(conn = conn, name = "2Fold", value = data, overwrite = TRUE)
                  
              }
              
              dbDisconnect(conn)
              return(paste0("Sites have been added to ", name))
              
          })
