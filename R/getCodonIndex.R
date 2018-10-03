#' Index genome according to provided exons
#'
#'
#' @description Author: Christopher Ward Index genome according to provided exons
#'
#' @param genome genome
#' @param exons GrangeList exons generated using \code{getFeatures} by passing \code{"gene:exons"} to \code{feature}
#'
#' @importFrom BSgenome getSeq
#'
#' @return An object of \code{CodonIndex}
#'
#'
#' @examples
#'
#' @export
#' @rdname getCodonIndex-methods


setGeneric("getCodonIndex", function(genome, exons){standardGeneric("getCodonIndex")})


#' @name getCodonIndex
#' @aliases getCodonIndex,CodonIndex-method
#' @rdname getCodonIndex-methods
#' @export
setMethod(f = "getCodonIndex", signature = "DNAStringSet",
          function(genome, exons, nCores){

            #initialize output
            out <- list(codonPositions = list(), RSCU = list())

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


              codList <- pbmclapply(seq(genes), mc.cores = nCores, function(x){

                iGene <- genes[[x]]
                refAnno <- exons[[x]]
                catGene <- unlist(iGene)

                if((length(catGene)/3)%%1 == 0 && !grepl("N", as.character(catGene), ignore.case = TRUE)){

                  # get codons cahnge to lower case
                  geneCodons <- tolower(codons(catGene))


                  # initialize aa sequence
                  aaCodons <-  geneCodons

                  residueRename <- lapply(1:length(lookUP), function(z){

                    aaCodons[grepl(lookUP[z], aaCodons)] <<- names(lookUP)[z]

                  })

                  rm(residueRename)


                  if(refAnno@strand@values == "+") {

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

                    seq(start(ex), end(ex))


                    })

                    positions <- do.call("c", positions)


                    fst <- positions[seq(1, length(positions), 3)]
                    snd <- positions[seq(2, length(positions), 3)]
                    trd <- positions[seq(3, length(positions), 3)]


                    data_frame(First = fst, Second = snd, Third = trd, Residue = aaCodons, Codon = geneCodons)

                  }

                }




              })


            names(codList) <- names(exons)
            codList <- Filter(Negate(is.null), codList)


            out[["codonPositions"]] <- codList


            # #calculate RSCU
            # codList <- lapply(1:length(catGenes), function(x){
            #
            #   #remove genes that do not have codons and check for missing data in the reference
            #   if((length(catGenes[[x]])/3)%%1 == 0 && !grepl("N", as.character(catGenes[[x]]), ignore.case = TRUE)){
            #     match <- t(sapply(FourFoldlookUP, grepl, dnaCodons))
            #
            #     redMatch <- Reduce(f = "|", split(match, rownames(match)))
            #
            #     dnaCodons4fold <-  dnaCodons[redMatch]
            #     count <- table(dnaCodons4fold)
            #
            #     RSCU <- lapply(FourFoldlookUP, function(z){
            #
            #       countSplit <- count[grepl(z, rownames(count))]
            #
            #       if(sum(countSplit) != 1) {
            #         RSCU <- lapply(countSplit, function(c){
            #
            #           (4*c)/sum(countSplit)
            #
            #         }) %>% unlist()
            #
            #       }
            #       else NULL
            #     }) %>% unlist()
            #
            #     df <- data_frame(codon = names(RSCU), RSCU = RSCU) %>% separate(codon, c("AA", "codon"))
            #
            #   }


            # })

            #set new args
            #args <- c(list(Class = "CodonIndex"))

            #do.call("new", args)

          })
