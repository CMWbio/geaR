#' build a codon sqlLite DB
#'
#'
#' @description 
#' Index genome according to codon positions in coding sequence. \cr
#' This can be useful when downstream analysis or output of alignements should only be applied to certain codons. \cr
#' For example calculating genetic distance on only 4-fold sites or outputting fasta files based on codon position.
#'
#' @param genome \code{character} or \code{DNAStringSet} \cr
#' Import fasta reference genome from file by passsing a \code{character} or provide a preloaded (using \code{Biostrings::readDNAStringSet()}) \code{DNAStringSet}
#' @param exons \code{GrangesList} \cr
#' Generated using \code{getFeatures()} by passing \code{"gene:cds"} to the  \code{feature} argument.
#' @param nCores \code{numeric}
#' Number of cores allcated.
#' @param position \code{character}
#' Options are \code{"all"}, default, or any combination of \code{c("first", "second", "third")}. \cr
#' Output will contain only specified codon positions.
#' @param fourFoldCodon \code{character}
#'  Options are \code{NULL}, \code{"only"} or \code{"remove"}. \cr
#' Determines what to do with 4-fold degenerate sites. \cr
#'  \cr
#' \code{NULL}, the default, will have no effect on the output, \code{"only"} will include only 4 fold degenerate sites in the output, \cr
#' \code{"exclude"} will output all sites (unless a position is specified) excluding 4-fold degenerate sites. \cr
#'  \cr
#' \code{"only"} is incompatable with specifying a position.
#'
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom future sequential
#' @importFrom Biostrings codons
#' @importFrom RSQLite SQLite
#' @importFrom RSQLite dbWriteTable
#' @importFrom RSQLite dbDisconnect
#' @importFrom RSQLite dbConnect
#'
#' @return An object of class \code{GRanges} containing codon ranges for positions specified.
#'
#'
#' @examples
#'
#' @export
#' @rdname buildCodonDB-methods


setGeneric("buildCodonDB", function(genome, exons, sqlDir = "./CodonStore.db", nCores = 1, fourFoldCodon = "include", ...){
  standardGeneric("buildCodonDB")
})

#' @aliases getCodonFeature,character
#' @export
setMethod("buildCodonDB", signature(genome = "character"),
          function(genome, exons, sqlDir = "./CodonStore.db", nCores = 1, fourFoldCodon = "include", ...){
            
            genome <- readDNAStringSet(genome)
            
            getCodonFeatures(genome, exons, nCores, sqlDir, fourFoldCodon, ...)
            
          })

#' @aliases getCodonFeature,DNAStringSet
#' @export
setMethod("buildCodonDB", signature = "DNAStringSet",
          function(genome, exons, sqlDir = "./CodonStore.db", nCores = 1, fourFoldCodon = "include"){
            
            # Residue lookup table
            lookUP <- c(Ala = "gc[tcag]",
                        Gly = "gg[tcag]",
                        Pro = "cc[tcga]",
                        Thr = "ac[tcga]",
                        Val = "gt[tagc]",
                        Arg4 = "cg[tagc]",
                        Leu4 = "ct[tagc]",
                        Ser4 = "tc[tagc]",
                        Arg2 = "ag[ag]",
                        Leu2 = "tt[ag]",
                        Ser2 = "ag[tc]",
                        Asn = "aa[tc]",
                        Asp = "ga[tc]",
                        Cys = "tg[tc]",
                        Gln = "ca[ag]",
                        Glu = "ga[ag]",
                        His = "ca[tc]",
                        Ile = "at[tca]",
                        Lys = "aa[ag]",
                        Met = "atg",
                        Phe = "tt[tc]",
                        Trp = "tgg",
                        Tyr = "ta[tc]",
                        Stp = "ta[ag]|tga"
            )
            
            
            #### add genes to exon grange this saves alot of memory
            #### why the fuck do DNAstringSets take up so much memory???
            exons <- .getGenes(genome, exons)
            
            
            codDF <- mclapply(seq(exons), mc.cores = nCores, .mapCodons, exons, lookUP)
            
            # codDF <- bind_rows(codDF)
            conn <- dbConnect(SQLite(), sqlDir)
            
            codDF2 <- map(seq(codDF)[1:10], function(x){
              d <- codDF[[x]]
              if(!is.null(d)){
                name <- d$gene[[1]]
                d$seqnames <- as.character(d$seqnames)
                d$strand <- as.character(d$strand)
                dbWriteTable(conn = conn, name = name, value = as.data.frame(d))
                name
              }
              
            })
            
            codDF2 <- unlist(codDF2)
            
            dbWriteTable(conn = conn, name = "genes", value = data.frame(genes = codDF2))
            # dbWriteTable(conn = conn, name = "First", value = as.data.frame(filter(codDF, codonPosition == "first")), overwrite = TRUE)
            # dbWriteTable(conn = conn, name = "Second", value = as.data.frame(filter(codDF, codonPosition == "second")), overwrite = TRUE)
            # dbWriteTable(conn = conn, name = "Third", value = as.data.frame(filter(codDF, codonPosition == "third")), overwrite = TRUE)
            codDF <- NULL
            dbDisconnect(conn)
            return(paste0("sqlDB built at ", sqlDir))
          })


.getGenes <- function(genome, exons){
  # get the genes from the genome using GRange list
  genes <- getSeq(genome, exons)
  exons@unlistData$seq <- as.character(unlist(genes))
  exons
}

.mapCodons <- function(x, exons, lookUP){
  
  refAnno <- exons[[x]]
  iGene <- DNAStringSet(refAnno$seq)
  
  catGene <- unlist(iGene)
  
  if((length(catGene)/3)%%1 == 0 & !grepl("N", as.character(catGene), ignore.case = TRUE)){
    
    # get codons cahnge to lower case
    geneCodons <- tolower(codons(catGene))
    
    # initialize aa sequence
    aaCodons <-  geneCodons
    
    aaCodons <- .replaceCodonNames(lookUP, aaCodons)
    
    positions <- .getCodonPositions(refAnno)
    
    fst <- positions[seq(1, length(positions), 3)]
    snd <- positions[seq(2, length(positions), 3)]
    trd <- positions[seq(3, length(positions), 3)]
    
    
    codDF <- .formatCodonStrands(refAnno, fst, snd, trd, geneCodons, aaCodons)
    
    # 
    # if(fourFoldCodon %in% c("only", "exclude")){
    #     
    #     fourFold <- c("Ala", "Gly", "Pro", "Thr", "Val", "Arg4", "Leu4", "Ser4")
    #     codDF$rownum <- 1:nrow(codDF)
    #     cod4Fold <- codDF[codDF$residue %in% fourFold,]
    #     
    #     
    # }
    # 
    # if(fourFoldCodon == "only"){
    #     codDF <- cod4Fold[colnames(cod4Fold) != "rownum"]
    #     
    # }
    # 
    # if(fourFoldCodon == "exclude"){
    #     cod4Fold <- cod4Fold[cod4Fold$codonPosition == "third",]
    #     
    #     codDF <- codDF[!codDF[["rownum"]] %in% cod4Fold[["rownum"]],]
    #     
    #     codDF <- codDF[colnames(codDF) != "rownum"]
    #     
    # }
    return(codDF)
  }
  
}

.formatCodonStrands <- function(refAnno, fst, snd, trd, geneCodons, aaCodons){
  if(refAnno@strand@values == "+") {
    codDF <- tibble(seqnames = refAnno@seqnames@values,start = fst, end = trd,
                    strand = refAnno@strand@values, gene = refAnno$gene[[1]],
                    residue = aaCodons, codon = geneCodons, first = fst, second = snd, third = trd, number = 1:length(aaCodons))
    codDF <- gather(codDF, "codonPosition", "start", c("first", "second", "third"))
    codDF["end"] <- codDF$start
    
    codDF <- codDF[order(codDF$start),]
  }
  if(refAnno@strand@values == "-") {
    codDF <- tibble(seqnames = refAnno@seqnames@values,
                    strand = refAnno@strand@values, gene = refAnno$gene[[1]],
                    residue = rev(aaCodons), codon = rev(geneCodons), first = fst, second = snd, third = trd, number = 1:length(aaCodons))
    
    codDF <- gather(codDF, "codonPosition", "start", c("first", "second", "third"))
    codDF["end"] <- codDF$start
    
    codDF <- codDF[order(-codDF$start),]
  }
  return(codDF)
  
}


.replaceCodonNames <- function(lookUP, aaCodons){
  residueRename <- lapply(1:length(lookUP), function(z){
    aaCodons[grepl(lookUP[z], aaCodons)] <<- names(lookUP)[z]
  })
  return(aaCodons)
}

.getCodonPositions <- function(refAnno){
  positions <- lapply(seq(refAnno), function(y){
    #get exon
    ex <- refAnno[y]
    
    ex <- seq(start(ex), end(ex))
    
    if(refAnno@strand@values == "-") ex <- rev(ex)
    
    ex
  })
  positions <- do.call("c", positions)
}

