#' Index genome according to provided exons
#'
#'
#' @description Author: Christopher Ward \cr
#' Converts GDS directly to fasta and outputs to a directory. \cr
#' Providing a range to the function will result in multiple files (number of ranges) being output with alignemnts within the ranges.
#'
#'
#' @param genome genome object of \code{DNAStringSet} class
#' @param exons GrangeList exons generated using \code{getFeatures} by passing \code{"gene:cds"} to \code{feature}
#' @param removeIndels removes indels, if false this will mess up the alignemtn in output fasta. concensus will also no longer work.1
#' @param fasta \code{DNAStringSet} the output of \code{Biostrings::readDNAStringSet()}. Input genome to extract features from, will align
#' automatically to the reference positions.
#'
#'
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings writeXStringSet
#' @importFrom Biostrings DNAStringSet
#' @importFrom msa msa
#'
#' @return  fasta in specified directory.
#'
#'
#' @examples
#'
#' @export
#' @rdname outputLociFasta


outputLociFasta <- function(GDS, loci, dir, pops, nCores = 1, ploidy = 2, alleles = "seperate", minSites = 0.1, removeIndels = TRUE, fasta = NULL){


  store <- mclapply(1:length(loci), mc.cores = nCores, function(locus){

   locus <- loci[[locus]]

    genoMat <- getGenotypes(GDS = GDS, pops = pops, locus = locus, minSites = minSites, nucleotide = TRUE, ploidy = ploidy, removeIndels = removeIndels)

    if(length(genoMat)){

    if(alleles == "concensus" & ploidy == 2){

      # look up table for diploids taken from BioStrings::IUPAC_CODE_MAP

      code <- c("W", "S", "M", "K", "R", "Y")

      codevec <- c("A", "C" ,"G" ,"T", rep(code, each = 2),"N")

      bases <- paste(c('A','C','G','T', 'A','T','C','G', 'A','C','G','T', 'A','G','C','T', "N"),
                     c('A','C','G','T', 'T','A','G','C', 'C','A','T','G', 'G','A','T','C', "N"),
                     sep = "/")

      names(codevec) <- bases

      names <- colnames(genoMat)
      pos <- rownames(genoMat)

      genoMat <- sapply(seq(from = 2, to = ncol(genoMat), by = 2), function(k){

        alleles <- paste(genoMat[,k-1], genoMat[,k], sep = "/")
        alleles <- codevec[as.vector(alleles)]

      })
      colnames(genoMat) <- unique(gsub("/.*", "", names))
      rownames(genoMat) <- pos

    }





    if("gene" %in% colnames(locus@elementMetadata)){

      st <- min(locus@ranges@start)
      ed <- min(end(locus))

      filename <- paste(as.character(GRanges(seqnames = locus@seqnames[1],
                                             strand = locus@strand[1],
                                             IRanges(start = st, end = ed))), collapse = ",")


      filename <- paste0(locus$gene[[1]], "_", filename, ".fasta")

    }

    else{

      filename <- paste(as.character(locus), collapse = ",")

      filename <- paste0(filename, ".fasta")

    }

      if(!is.null(fasta)) {
        ref <- getSeq(fasta, locus)
        ref <- unlist(ref)
        ref <- DNAStringSet(ref)
        ref@ranges@NAMES <- "Ref"
        ref <- as_tibble(genos[[1]])
        ref <- rownames_to_column(ref, "position")

        colnames(ref) <- c("position", "Ref")

        genoPos <- rownames(genoMat)
        genoMat <- as.data.frame(genoMat)
        rownames(genoMat) <- genoPos
        genoMat <- rownames_to_column(genoMat, "position")

        genoMat <- left_join(ref, genoMat, by = "position")
        genoMat <- select(genoMat, -position)
        genoMat[is.na(genoMat)] <- "N"

        genoMat <- as.matrix(genoMat)

      }

    genoMat <- t(genoMat)
    genos <- split(genoMat, factor(rownames(genoMat), levels = rownames(genoMat)))
    genos <- sapply(genos, paste, collapse="")
    genos <- DNAStringSet(genos)

    if(all(locus@strand == "-")) {

      genos <- reverseComplement(genos)
    }
    if(align) writeXStringSet(genos@unmasked, paste0(dir, "/", filename))
    else  writeXStringSet(genos, paste0(dir, "/", filename))



 }

  })


  }
