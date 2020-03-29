#' getPositionFromDB
#'
#' @description extracts a partition scheme from a codon DB 
#'
#' @details 
#' Accepts both GRange and SQL codon databases generated using \code{buildCodonDB}
#' 
#' @param DB \code{GRange or character}. Either the path to a SQL database or GRange object 
#' containing codon information.
#' @param positions \code{character} vector containing any of code{c("first", 
#' "second", "third")}
#'
#' @return a GRange object containing specified codon positions.
#'
#'
#' @rdname getPositionFromDB-methods
#' @export
setGeneric("getPositionFromDB",function(DB, positions){
    standardGeneric("getPositionFromDB")
})

#' @aliases getPositionFromDB,character
#' @export
setMethod("getPositionFromDB", signature(DB = "character"),
          function(DB, positions){
              
              conn <- dbConnect(SQLite(), DB)
              genes <- dbGetQuery(conn, paste0("SELECT * FROM genes"))[[1]]
              
              
              data <- mclapply(genes, mc.cores = nCores, function(x){
                  # x <- genes[[x]]
                  
                  gene <- dbGetQuery(conn, paste0("SELECT * FROM '", x, "'"))
                  
                  gene <- gene[gene$codonPosition %in% positions,]
                  
                  gene <- gene[colnames(gene) %in% c("seqnames", "start", "end", "strand", "gene")]
                  
              })
              
              data <- bind_rows(data)
              
              data$Name <-  data$gene
              
              makeGRangesListFromDataFrame(data, split.field = "Name", keep.extra.columns=TRUE)
          })

#' @aliases getPositionFromDB,character
#' @export
setMethod("getPositionFromDB", signature(DB = "GRangesList"),
          function(DB, positions){
              
              data <- mclapply(seq(DB), mc.cores = nCores, function(x){
                  gene <- DB[[x]]
                  gene <- as_tibble(gene)
                  gene <- gene[gene$codonPosition %in% positions,]
                  
                  gene <- gene[colnames(gene) %in% c("seqnames", "start", "end", "strand", "gene")]
              })
              
              data <- bind_rows(data)
              
              data$Name <-  data$gene
              
              makeGRangesListFromDataFrame(data, split.field = "Name", keep.extra.columns=TRUE)
              
          })