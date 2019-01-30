#' Index genome according to provided exons
#'
#'
#' @description Author: Christopher Ward \cr
#' Index genome according to codon positions in coding sequence. \cr
#' This can be useful when downstream analysis or output of alignements should only be applied to certain codons. \cr
#' For example calculating genetic Tajima's D on only 4-fold sites or outputting fasta files based on codon position.
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
#'
#' @return An object of class \code{GRanges} containing codon ranges for positions specified.
#'
#'
#' @examples
#'
#' @export
#' @rdname getCodonFeatures-methods


setGeneric("getCodonFeatures", function(genome, exons, nCores = 1, position = "all",  fourFoldCodon = "include", ...){
  standardGeneric("getCodonFeatures")
})

#' @aliases getCodonFeature,character
#' @export
setMethod("getCodonFeatures", signature(genome = "character"),
          function(genome, exons, nCores = 1, position = "all",  fourFoldCodon = "include", ...){

            genome <- readDNAStringSet(genome)

            getCodonFeatures(genome, exons, nCores, position,  fourFoldCodon, ...)

          })

#' @aliases getCodonFeature,DNAStringSet
#' @export
setMethod("getCodonFeatures", signature = "DNAStringSet",
          function(genome, exons, nCores = 1, position = "all",  fourFoldCodon = "include"){

            if(position == "all") position <- c("first", "second", "third")

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




            # get the genes from the genome using GRange list
            genes <- BSgenome::getSeq(genome, exons)


            codList <- mclapply(seq(exons), mc.cores = nCores, function(x){

              iGene <- genes[[x]]
              refAnno <- exons[[x]]

              catGene <- unlist(iGene)

              if((length(catGene)/3)%%1 == 0 & !grepl("N", as.character(catGene), ignore.case = TRUE)){

                # get codons cahnge to lower case
                geneCodons <- tolower(codons(catGene))

                # initialize aa sequence
                aaCodons <-  geneCodons

                residueRename <- lapply(1:length(lookUP), function(z){

                  aaCodons[grepl(lookUP[z], aaCodons)] <<- names(lookUP)[z]

                })

                rm(residueRename)


                #if(refAnno@strand@values == "+") {

                # get teh exon widths for the gene
                #exonWidths <- iGene@ranges@width / 3

                # sum to get the sfinal position in each exon codon position
                #lastCodon <- cumsum(exonWidths)
                # get the first codon numbers
                #firstCodon <- as.integer(lastCodon - exonWidths + 1)

                #lastCodon <- as.integer(lastCodon)

                positions <- lapply(seq(refAnno), function(y){

                  #get exon
                  ex <- refAnno[y]

                  ex <- seq(start(ex), end(ex))

                  if(refAnno@strand@values == "-") ex <- rev(ex)

                  ex


                })

                positions <- do.call("c", positions)


                fst <- positions[seq(1, length(positions), 3)]
                snd <- positions[seq(2, length(positions), 3)]
                trd <- positions[seq(3, length(positions), 3)]


                if(refAnno@strand@values == "+") {
                  codDF <- tibble(seqnames = refAnno@seqnames@values,start = fst, end = trd,
                                      strand = refAnno@strand@values, gene = exons[[x]]$Name[[1]],
                                      residue = aaCodons, codon = geneCodons, first = fst, second = snd, third = trd)
                  codDF <- gather(codDF, "codonPosition", "start", c("first", "second", "third"))
                  codDF["end"] <- codDF$start

                  codDF <- codDF[order(codDF$start),]
                }
                if(refAnno@strand@values == "-") {
                  codDF <- tibble(seqnames = refAnno@seqnames@values,
                                      strand = refAnno@strand@values, gene = exons[[x]]$Name[[1]],
                                      residue = rev(aaCodons), codon = rev(geneCodons), first = fst, second = snd, third = trd)

                  codDF <- gather(codDF, "codonPosition", "start", c("first", "second", "third"))
                  codDF["end"] <- codDF$start

                  codDF <- codDF[order(-codDF$start),]
                }


                if(fourFoldCodon %in% c("only", "exclude")){

                  fourFold <- c("Ala", "Gly", "Pro", "Thr", "Val", "Arg4", "Leu4", "Ser4")
                  codDF$rownum <- 1:nrow(codDF)
                  cod4Fold <- codDF[codDF$residue %in% fourFold,]


                }

                if(fourFoldCodon == "only"){
                  position <- "third"

                  codDF <- cod4Fold[colnames(cod4Fold) != "rownum"]

                }

                if(fourFoldCodon == "exclude"){
                  cod4Fold <- cod4Fold[cod4Fold$codonPosition == "third",]

                  codDF <- codDF[!codDF[["rownum"]] %in% cod4Fold[["rownum"]],]

                  codDF <- codDF[colnames(codDF) != "rownum"]

                }




                codDF <- codDF[codDF$codonPosition %in% position,]

                codGR <- makeGRangesFromDataFrame(codDF, keep.extra.columns = TRUE)
              }


            })

            codList <- Filter(Negate(is.null), codList)
            GRangesList(codList)

          })
