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
#' @importFrom RSQLite dbGetQuery
#'
#' @examples
#'
#' @export
#' @rdname valicate4FoldCodons

validate4FoldCodons <- function(GDS, sqlDir, nCores, pops){
    
    fourFold <- c("Ala", "Gly", "Pro", "Thr", "Val", "Arg4", "Leu4", "Ser4")
    
    conn <- dbConnect(SQLite(), sqlDir)
    genes <- dbGetQuery(conn, paste0("SELECT * FROM genes"))[[1]]

    
    samples <- pops$Sample
    seqSetFilter(object = GDS, sample.id = samples)
    
    data <- mclapply(genes, mc.cores = nCores, function(x){
        # x <- genes[[x]]

        gene <- dbGetQuery(conn, paste0("SELECT * FROM '", x, "'"))
 
        .validAF(gene, fourFold, pops, GDS)
 
    })

    data <- bind_rows(data)
    dbWriteTable(conn = conn, name = "4Degen", value = data, overwrite = TRUE)
    dbDisconnect(conn)
    return(paste0("Sites have been added to ", sqlDir))
}

.validAF <- function(gene, fourFold, pops, GDS){
    
    F_S <- gene[gene$residue %in% fourFold & gene$codonPosition %in% c("first", "second"),]
    F_S <- makeGRangesFromDataFrame(F_S, keep.extra.columns = TRUE)

    Td <- gene[gene$residue %in% fourFold & gene$codonPosition %in% c("third"),]
    
    seqSetFilter(object = GDS, variant.sel = F_S)
    
    pos <- seqGetData(GDS, var.name = "position")
    alt <- seqAlleleFreq(GDS)
    pos <- pos[alt == 1]
    
    trueF_S <- F_S[F_S@ranges@start %in% pos,] 
    trueF_S <- as_tibble(trueF_S)
    trueF_S <- group_by(trueF_S, number)
    trueF_S <- filter(trueF_S, length(number) == 2) 
    trueTd <- filter(Td, number %in% trueF_S$number)
    return(trueTd)
}

