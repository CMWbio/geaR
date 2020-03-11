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
#' @importFrom purrr map
#' @importFrom furrr future_map
#' @importFrom future plan
#'
#' @return A \code{GRangesList} object where each element is a genomic feature
#'
#'
#' @examples
#'
#' @name makeFeatures
#' @rdname makeFeatures-methods
#' @export

setGeneric("makeFeatures", function(gffName,  feature, nCores = 1,
                                   longestIsoform, includeRange, geneIDField, ...){
  standardGeneric("makeFeatures")
})


#' @aliases makeFeatures,which-Granges
#' @export
setMethod("makeFeatures", signature(gffName = "character",
                                   includeRange = "GRanges"),
          function(gffName, feature, nCores,
                   longestIsoform, includeRange, geneIDField){
            G_range <- import.gff(gffName, which = includeRange)
            makeFeatures(gffName = G_range, feature = feature, nCores = nCores, longestIsoform = longestIsoform, geneIDField = geneIDField)
          })

#' @aliases makeFeatures,which-GrangesList
#' @export
setMethod("makeFeatures", signature(gffName = "character",
                                   includeRange = "GRangesList"),
          function(gffName, feature, nCores,
                   longestIsoform, includeRange, geneIDField){
            G_range <- import.gff(gffName, which = includeRange)
            makeFeatures(gffName = G_range, feature = feature, nCores = nCores, longestIsoform = longestIsoform, geneIDField = geneIDField)
          })

#' @aliases makeFeatures,which
#' @export
setMethod("makeFeatures", signature(gffName = "character"),
          function(gffName, feature, nCores,
                   longestIsoform, geneIDField){
            G_range <- import.gff(gffName)
            makeFeatures(gffName = G_range, feature = feature, nCores = nCores, longestIsoform = longestIsoform, geneIDField = geneIDField)
          })



#' @aliases makeFeatures,import
#' @export
setMethod("makeFeatures", signature = "GRanges",
          function(gffName, feature, nCores,
                   longestIsoform, geneIDField){

            # check ig gene is specified in the feature parameter
            if(grepl("gene", feature)){
              
              if(missing(geneIDField)) geneIDField <- "ID"
              gffName$gene <- gffName@elementMetadata[[geneIDField]]
              gffName@elementMetadata <- gffName@elementMetadata[,colnames(gffName@elementMetadata) %in% c("type", "ID", "Name", "gene", "Parent")]


              # filter all GRanges to contain only those with type == gene
              genes <- gffName[gffName$type == "gene",]
              mRNA <- gffName[gffName$type == "mRNA"]

              # get ID from Granges using the ID field
              geneID <- genes$gene

              # make into GRangesList containing only genes
              genesList <- GRangesList(
                split(genes,
                      geneID))

              # if only genes are specified then return the genes GRangeList
              if(feature == "gene") {
                return(genesList)
              }
              # Continue if you also want exons
              if(feature == "gene:exons"){

                # get all exons from all GRanges
                all <- gffName[gffName$type == "exon",]
                all@elementMetadata <- all@elementMetadata[colnames(all@elementMetadata) %in% c("gene", "Parent")]

                if(length(unlist(mRNA$Parent))){
                  mRNA_df <- tibble::tibble(geneID = unlist(mRNA$Parent), mrnaID = mRNA$ID)

                  gene_df <- tibble::tibble(geneName = genes$gene, geneID = genes$ID)

                  comb_df <- dplyr::full_join(mRNA_df, gene_df)

                  comb_df$mrnaID[is.na(comb_df$mrnaID)] <- comb_df$geneID[is.na(comb_df$mrnaID)]

                  mRNAlookup <- comb_df$geneID
                  names(mRNAlookup) <- comb_df$mrnaID

                  gNames <- unlist(all$Parent)
                  gNames <- unname(mRNAlookup[gNames])

                  geneLookup <-  comb_df$geneName
                  names(geneLookup) <- comb_df$geneID

                  gNames <- unname(geneLookup[gNames])

                  all$gene <- gNames
                } else all$gene <- unlist(all$Parent)

                all <- split(all, all$gene)

                # get mRNA using the gene ID field


                ExonName <- names(all)

                plan(multiprocess, workers = nCores)
                all <- future_map(seq_along(all), .f = .mapFun, all, longestIsoform)

                all <- dplyr::bind_rows(all)
                all <- makeGRangesListFromDataFrame(all, split.field = "gene", keep.extra.columns = TRUE)            
                return(all)
              }
              if(feature == "gene:cds"){

                # get all exons from all GRanges
                all <- gffName[gffName$type == "CDS",]

                all@elementMetadata <- all@elementMetadata[colnames(all@elementMetadata) %in% c("gene", "Parent")]

                #make lookup table for CDS parent tracking

                if(length(unlist(mRNA$Parent))){
                  mRNA_df <- tibble::tibble(geneID = unlist(mRNA$Parent), mrnaID = mRNA$ID)

                  gene_df <- tibble::tibble(geneName = genes$gene, geneID = genes$ID)

                  comb_df <- dplyr::full_join(mRNA_df, gene_df)

                  comb_df$mrnaID[is.na(comb_df$mrnaID)] <- comb_df$geneID[is.na(comb_df$mrnaID)]

                  mRNAlookup <- comb_df$geneID
                  names(mRNAlookup) <- comb_df$mrnaID

                  gNames <- unlist(all$Parent)
                  gNames <- unname(mRNAlookup[gNames])

                  geneLookup <-  comb_df$geneName
                  names(geneLookup) <- comb_df$geneID

                  gNames <- unname(geneLookup[gNames])

                  all$gene <- gNames
                } else all$gene <- unlist(all$Parent)


                all <- split(all, all$gene)

                plan(multiprocess, workers = nCores)

                all <- future_map(seq_along(all), .f = .mapFun, all, longestIsoform)

                all <- dplyr::bind_rows(all)
                all <- makeGRangesListFromDataFrame(all, split.field = "gene", keep.extra.columns = TRUE)
                

                #CDS <- Filter(length, CDS)

                return(all)
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

.mapFun <- function(y, all, longestIsoform){
    
    rec <- all[[y]]
    
    mRNA <- split(rec, as.character(rec$Parent))
    
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
        longest <- which(isoformLengths == max(isoformLengths))[1]
        mRNA <- mRNA[[longest]]} else mRNA <- mRNA[[1]]
    
    tibble::as_tibble(mRNA)
    
}


