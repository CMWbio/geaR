#' Get loci from gff
#'
#' @description Creates feature loci
#'
#' @details Authour: Chris Ward
#'
#' @param gffName path to gff to extract features from
#' @param contigMD data frame with contig metadata. Contains two columns: \code{ID} with scaffold names, \code{length} with scaffold length
#' @param feature \code{list} or \code{GRangesList}. Loci to import genotypes for. c("gene", "gene:exon", "gene:cds", "pseudogene", "lncRNA", "intergenic")
#' @param minSites \code{numeric} minimum number of sites as a proportion of loci length. Default 0.5 (ie 50 percent)
#' @param nCores \code{numeric} number of cores to run in parallel
#' @param longestIsoform \code{logical} only effects feature = "gene:exon". by default select the first entry for each gene. If \code{TRUE} will select the longest isoform for that gene. If there is a biotype field in the GFF it will only deal with those labeled as 'protein coding.' Selecting the longest isoform may cause be problematic if you do not have the biotype.
#'
#' @importFrom rtracklayer import.gff
#' @import pbmcapply
#'
#' @return A \code{list} of gene GRanges
#'
#'
#' @examples
#'
#' @export
#' @rdname getFeatures

getFeatures <- function(gffName, contigMD, feature = "gene", nCores = 1, longestIsoform = FALSE){


  # read in the gff3
  allGR <- import.gff(gffName)

  # check ig gene is specified in the feature parameter
  if(grepl("gene", feature)){

    # filter all GRanges to contain only those with type == gene
    genes <- allGR[allGR$type == "gene",]
    # make into GRangesList containing only genes
    genes <- GRangesList(split(genes, genes$Name))

    # if only genes are specified then return the genes GRangeList
    if(feature == "gene") return(genes)

    # Continue if you also want exons
    if(feature == "gene:exons"){

      # get all exons from all GRanges
      allExons <- allGR[allGR$type == "exon",]

      # Start making Coding sequences, will select the longest isoform or the first entry for each gene
      exons <- GRangesList(mclapply(seq_along(genes), mc.cores = nCores, function(x){

        #get Grange using index position
        gr <- genes[[x]]
        # get mRNA using the gene ID field
        mRNA <- subset(allGR, allGR$Parent == gr$ID)

        # select longest isoform
        if(longestIsoform){
          #check for biotype to only select isoforms with 'protien_coding' and 'TEC'
          if("biotype" %in% names(mRNA@elementMetadata)){

            if(length(mRNA$biotype) != 1) mRNA <- subset(mRNA, biotype %in% c("protein_coding", "TEC"))

          }

          # get the max length sequence
          mRNA <- mRNA[which(mRNA@ranges@width == max(mRNA@ranges@width))]

        }else{

          mRNA <- mRNA[1]

        }

        # subset exons by those within selected mRNA
        subset(allExons, allExons$Parent == mRNA$ID)

      }))

      #name the GRanges List containing coding sequence
      names(exons) <- names(genes)
      return(exons)
    }
    if(feature == "gene:cds"){

      # get all exons from all GRanges
      allCDS <- allGR[allGR$type == "CDS",]

      # Start making Coding sequences, will select the longest isoform or the first entry for each gene
      CDS <- GRangesList(mclapply(seq_along(genes), mc.cores = nCores, function(x){

        #get Grange using index position
        gr <- genes[[x]]
        # get mRNA using the gene ID field
        mRNA <- subset(allGR, allGR$Parent == gr$ID)

         #select longest isoform
         if(longestIsoform){
           #check for biotype to only select isoforms with 'protien_coding' and 'TEC'
           #if("biotype" %in% names(mRNA@elementMetadata)){

        #     if(length(mRNA$biotype) != 1) mRNA <- subset(mRNA, biotype %in% c("protein_coding", "TEC"))
        #
        #   }
        #
           # get the max length sequence
           mRNA <- mRNA[which(mRNA@ranges@width == max(mRNA@ranges@width))]

           if(length(unique(mRNA$ID)) > 1) mRNA <- mRNA[1]

         }else{

          mRNA <- mRNA[1]

         }

        # subset exons by those within selected mRNA
      subset(allCDS, allCDS$Parent == mRNA$ID)

      }))

      #name the GRanges List containing coding sequence
      names(CDS) <- names(genes)
      CDS <- Filter(length, CDS)

      return(CDS)
    }


  }



  #get the psuedogenes out
  if(feature == "pseudogene"){

    psGenes <- allGR[allGR$type == "pseudogene"]
    psGenes <- GRangesList(split(psGenes, psGenes$Name))
    return(psGenes)

  }

  #get the lncRNA out
  if(feature == "lncRNA"){

    lncRNA <- allGR[allGR$type == "lnc_RNA"]
    lncRNA <- GRangesList(split(lncRNA, lncRNA$Name))
    return(lncRNA)

  }

  # get the intergenic sequences out
  if(feature == "intergenic"){

    #get contig range from the contigMD object
    contigRange <- makeGRangesFromDataFrame(data_frame(chr = contigMD[["ID"]], start = 1, strand = ".", end = contigMD[["length"]]))
    # remove any chromosome annotations
    noChr <- allGR[allGR$type != "chromosome"]
    # join all annotations that overlap
    noChr <- disjoin(noChr)
    # get hits on each sequence and get overlap
    hits <- findOverlaps(contigRange, noChr)
    overlapHits <- extractList(noChr, as(hits, "List"))
    #use hits to get intergenic regions
    intergenic <- psetdiff(contigRange, grl)

    return(intergenic)

  }






}



