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
#' @param removeIndels removes indels, if false this will mess up the alignemtn in output fasta. concensus will also no longer work.
#' @param
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


outputLociFasta <- function(GDS, loci, dir, pops, nCores = 1, ploidy = 2, alleles = "seperate", minSites = 0.1, removeIndels = TRUE, fasta = NULL, align = FALSE){


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

      genoMat <- sapply(seq(from = 2, to = ncol(genoMat), by = 2), function(k){

        alleles <- paste(genoMat[,k-1], genoMat[,k], sep = "/")
        alleles <- codevec[as.vector(alleles)]

      })
      colnames(genoMat) <- unique(gsub("/.*", "", names))
      rownames(genoMat) <- 1:nrow(genoMat)

    }

    genoMat <- t(genoMat)

    genos <- split(genoMat, rownames(genoMat))

    if("Name" %in% colnames(locus@elementMetadata)){

      filename <- paste0(locus$Name[[1]], ".fasta")

    }

    else{

      filename <- paste(as.character(locus), collapse = ",")
      filename <- paste0(filename, ".fasta")

      }


    genos <- sapply(genos, paste, collapse="")
    genos <- DNAStringSet(genos)

    if(all(locus@strand == "-")) {

      genos <- reverseComplement(genos)
    }

    if(!is.null(fasta)) {
      ref <- getSeq(fasta, locus)
      ref <- unlist(ref)
      ref <- DNAStringSet(ref)
      ref@ranges@NAMES <- "Ref"
      genos <- c(ref, genos)

      if(align){

      genos <- msa(genos, order = "input")
      }
    }


    writeXStringSet(genos@unmasked, paste0(dir, "/", filename))

 }

  })


  }
