#' Validate 4fold degenerate codons
#'
#'
#' @description This function will validate 4 fold degenerate codons using a provided SNP dataset
#' You only need to run this if you are interested in 4 fold degenerate sights 
#' Usually used for their quasi neutrality
#'
#'
#' @return An object of class \code{GRanges} containing codon ranges for positions specified.
#'
#' @param GDS an open GDS file with genptypes
#' @param input \code{character} or \code{GRangesList}. Path to sql database or GRange list containing codon index
#' @param nCores \code{numeric} number of cores
#' @param pops \code{data.frame} population data frame
#'
#' @importFrom RSQLite dbGetQuery
#'
#' @examples
#'
#' @export
#' @rdname validate4FoldCodons


setGeneric("validate4FoldCodons", function(GDS, input, nCores, pops){
    standardGeneric("validate4FoldCodons")
})

#' @aliases validate4FoldCodons,character
#' @export
setMethod("validate4FoldCodons", signature(input = "character"),
          function(GDS, input, nCores, pops){
              
              fourFold <- c("Ala", "Gly", "Pro", "Thr", "Val", "Arg4", "Leu4", "Ser4")
              
              conn <- dbConnect(SQLite(), input)
              genes <- dbGetQuery(conn, paste0("SELECT * FROM genes"))[[1]]
              
              
              samples <- pops$Sample
              seqSetFilter(object = GDS, sample.id = samples)
              
              data <- mclapply(genes, mc.cores = nCores, function(x){
                  # x <- genes[[x]]
                  
                  gene <- dbGetQuery(conn, paste0("SELECT * FROM '", x, "'"))
                  
                  if(length(gene)) df <- .validAF(gene, fourFold, pops, GDS)
                  else NULL
                  if(length(df)) df
                  else NULL
                  
              })
              
              data <- dplyr::bind_rows(data)
              
              dbWriteTable(conn = conn, name = "4Degen", value = data, overwrite = TRUE)
              dbDisconnect(conn)
              return(paste0("Sites have been added to ", input))
          })

#' @aliases validate4FoldCodons,character
#' @export
setMethod("validate4FoldCodons", signature(input = "GRangesList"),
          function(GDS, input, nCores, pops){
              
              fourFold <- c("Ala", "Gly", "Pro", "Thr", "Val", "Arg4", "Leu4", "Ser4")
              
              samples <- pops$Sample
              seqSetFilter(object = GDS, sample.id = samples)
              
              data <- mclapply(seq(input), mc.cores = nCores, function(x){
                   gene <- input[[x]]
                  
                  if(length(gene)) df <- .validAF(gene, fourFold, pops, GDS)
                  else NULL
                  if(length(df)) df
                  else NULL
                  
              })
              
              data <- dplyr::bind_rows(data)
              
              data$Name <-  data$gene
              
              makeGRangesListFromDataFrame(data, split.field = "Name", keep.extra.columns=TRUE)

          })




.validAF <- function(gene, fourFold, pops, GDS){
    
    
    if(class(gene) == "GRanges"){
        
        F_S <- gene[gene$residue %in% fourFold & gene$codonPosition %in% c("first", "second"),]
        Td <- as.data.frame(gene[gene$residue %in% fourFold & gene$codonPosition %in% c("third"),])
    }
    else{
        F_S <- gene[gene$residue %in% fourFold & gene$codonPosition %in% c("first", "second"),]
        F_S <- makeGRangesFromDataFrame(F_S, keep.extra.columns = TRUE)
        
        Td <- gene[gene$residue %in% fourFold & gene$codonPosition %in% c("third"),]
    }
    
        seqSetFilter(object = GDS, variant.sel = F_S)
        
        pos <- seqGetData(GDS, var.name = "position")
        if(length(pos)){
        alt <- seqAlleleFreq(GDS)
        pos <- pos[alt == 1]
        
        trueF_S <- F_S[F_S@ranges@start %in% pos,] 
        trueF_S <- as_tibble(trueF_S)
        trueF_S <- group_by(trueF_S, number)
        trueF_S <- filter(trueF_S, length(number) == 2) 
        trueTd <- filter(Td, number %in% trueF_S$number)
        return(trueTd)
    } else NULL
    
}

