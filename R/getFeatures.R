#' getFeatures
#'
#' @description Creates feature loci from a gff file
#'
#' @details Authour: Chris Ward \cr
#' Get genomic features such as 'genes', 'coding sequence', or 'exons' from gff.
#'
#' @param gffName \code{character} \cr
#' Path to gff to extract features from
#' @param contigMD \code{data.frame} \cr
#' Contig metadata. Contains two columns: \code{ID} with scaffold names, \code{length} with scaffold length. \cr
#' Used to select intergenic regions.
#' @param feature \code{character}. \cr
#' Loci to import genotypes for. c("gene", "gene:exon", "gene:cds", "pseudogene", "lncRNA", "intergenic")
#' @param nCores \code{numeric} \cr
#' Number of cores to run in parallel
#' @param includeRange \code{GRanges} or \code{GRangesList} \cr
#' Target ranges to extract features over.
#' @param longestIsoform \code{logical} \cr
#' Only effects feature = "gene:exon". \cr
#' By default function will select the first entry for each gene. \cr
#' If \code{TRUE} is passed, will select the longest isoform for that gene. \cr
#' If there is a \code{"biotype"} field in the GFF it will only deal with those labeled as 'protein coding.' \cr
#' Selecting the longest isoform may be problematic if the \code{"biotype"} is missing.
#' @param geneIdField \code{character} \cr
#' Name of the gene id field in the gff, \code{"Name"} by default.
#'
#' @import rtracklayer
#' @import pbmcapply
#'
#' @return A \code{GRangesList} object where each element is a genomic feature
#'
#'
#' @examples
#'
#' @name getFeatures
#' @rdname getFeatures-methods
#' @export

setGeneric("getFeatures", function(gffName,  feature = "gene:cds", nCores = 1,
                                   longestIsoform = FALSE, includeRange, geneIdField = "Name", ...){
  standardGeneric("getFeatures")
})


#' @aliases getFeatures,which-Granges
#' @export
setMethod("getFeatures", signature(gffName = "character",
                                   includeRange = "GRanges"),
          function(gffName, feature = "gene:cds", nCores = 1, longestIsoform = FALSE, includeRange, ...){
            G_range <- import.gff(gffName, which = includeRange)
            getFeatures(gffName = G_range, feature, nCores, longestIsoform, ...)
          })

#' @aliases getFeatures,which-GrangesList
#' @export
setMethod("getFeatures", signature(gffName = "character",
                                   includeRange = "GRangesList"),
          function(gffName, feature = "gene:cds", nCores = 1, longestIsoform = FALSE, includeRange, ...){
            G_range <- import.gff(gffName, which = includeRange)
            getFeatures(gffName = G_range, feature, nCores, longestIsoform, ...)
          })

#' @aliases getFeatures,which
#' @export
setMethod("getFeatures", signature(gffName = "character"),
          function(gffName, feature = "gene:cds", nCores = 1, longestIsoform = FALSE, ...){
            G_range <- import.gff(gffName)
            getFeatures(gffName = G_range, feature, nCores, longestIsoform, ...)
          })



#' @aliases getFeatures,import
#' @export
setMethod("getFeatures", signature = "GRanges",
          function(gffName, feature = "gene:cds", nCores = 1, longestIsoform = FALSE, includeRange, ...){

            # check ig gene is specified in the feature parameter
            if(grepl("gene", feature)){


              # filter all GRanges to contain only those with type == gene
              genes <- gffName[gffName$type == "gene",]

              # get ID from Granges using the ID field
              geneID <- genes@elementMetadata[colnames(genes@elementMetadata) == geneIdField][[1]]

              # make into GRangesList containing only genes
              genes <- GRangesList(
                split(genes,
                      geneID))

              # if only genes are specified then return the genes GRangeList
              if(feature == "gene") return(genes)

              # Continue if you also want exons
              if(feature == "gene:exons"){

                # get all exons from all GRanges
                allExons <- gffName[gffName$type == "exon",]

                # Start making Coding sequences, will select the longest isoform or the first entry for each gene
                exons <- GRangesList(mclapply(seq_along(genes), mc.cores = nCores, function(x){

                  #get Grange using index position
                  gr <- genes[[x]]
                  geneName <- names(genes)[x]

                  # get mRNA using the gene ID field
                  mRNA <- subset(gffName, gffName$Parent == gr$ID)

                  # select longest isoform
                  if(longestIsoform){
                    #check for biotype to only select isoforms with 'protien_coding' and 'TEC'
                    if("biotype" %in% names(mRNA@elementMetadata)){

                      if(length(mRNA$biotype) != 1) mRNA <- subset(mRNA, biotype %in% c("protein_coding", "TEC"))

                    }

                    # get the max length sequence
                    mRNA <- mRNA[which(mRNA@ranges@width == max(mRNA@ranges@width))]


                  }

                  mRNA <- mRNA[1]



                  # subset exons by those within selected mRNA
                  subset(allExons, allExons$Parent == mRNA$ID)

                }))

                #name the GRanges List containing coding sequence
                names(exons) <- names(genes)
                return(exons)
              }
              if(feature == "gene:cds"){

                # get all exons from all GRanges
                allCDS <- gffName[gffName$type == "CDS",]

                cdsList <- GRangesList(split(allCDS, unlist(allCDS$Parent)))


                # Start making Coding sequences, will select the longest isoform or the first entry for each gene
                CDS <- GRangesList(mclapply(seq_along(cdsList), mc.cores = nCores, function(x){



                  #get Grange using index position
                  gr <- cdsList[[x]]

                  noNA <- gffName[!is.na(gffName$ID),]
                  transcript <- noNA[noNA$ID ==  gr$Parent[[1]],]
                  gene <- noNA[noNA$ID == transcript$Parent[[1]],]

                  gr$Name <- gene$Name

                  gr


                }))

                # get mRNA using the gene ID field
                CDS <- split(CDS@unlistData, as.character(CDS@unlistData$Name))



                CDS <- GRangesList(mclapply(seq_along(CDS), mc.cores = nCores, function(y){

                  rec <- CDS[[y]]

                  mRNA <- split(rec, rec$ID)

                  #select longest isoform
                  if(longestIsoform){
                    #check for biotype to only select isoforms with 'protien_coding' and 'TEC'
                    #if("biotype" %in% names(mRNA@elementMetadata)){

                    #     if(length(mRNA$biotype) != 1) mRNA <- subset(mRNA, biotype %in% c("protein_coding", "TEC"))
                    #
                    #   }
                    #

                    isoformLengths <- unlist(lapply(mRNA, function(x){

                      sum(width(x))
                    }))

                    # get the max length sequence, extract first element in longest if there are multiple longest
                    longest <- which(isoformLengths == max(isoformLengths))
                    mRNA <- mRNA[[longest]]}
                    mRNA <- mRNA[[1]]


                }))



                #name the GRanges List containing coding sequence
                CDS <- Filter(length, CDS)

                list_elt_seqnames <- as.character(runValue(seqnames(CDS)))
                list_elt_seqnames <- factor(list_elt_seqnames, levels=seqlevels(CDS))
                list_elt_smallest_start <- min(start(CDS))
                oo <- order(as.integer(list_elt_seqnames), list_elt_smallest_start)
                CDS <- CDS[oo]

                return(CDS)
              }


            }



            #get the psuedogenes out
            if(feature == "pseudogene"){

              psGenes <- gffName[gffName$type == "pseudogene"]
              psGenes <- GRangesList(split(psGenes, psGenes$Name))
              return(psGenes)

            }

            #get the lncRNA out
            if(feature == "lncRNA"){

              lncRNA <- gffName[gffName$type == "lnc_RNA"]
              lncRNA <- GRangesList(split(lncRNA, lncRNA$Name))
              return(lncRNA)

            }

            # get the intergenic sequences out
            if(feature == "intergenic"){

              #get contig range from the contigMD object
              contigRange <- makeGRangesFromDataFrame(data_frame(chr = contigMD[["ID"]], start = 1, strand = ".", end = contigMD[["length"]]))
              # remove any chromosome annotations
              noChr <- gffName[gffName$type != "chromosome"]
              # join all annotations that overlap
              noChr <- disjoin(noChr)
              # get hits on each sequence and get overlap
              hits <- findOverlaps(contigRange, noChr)
              overlapHits <- extractList(noChr, as(hits, "List"))
              #use hits to get intergenic regions
              intergenic <- psetdiff(contigRange, overlapHits)

              return(intergenic)

            }



          })
