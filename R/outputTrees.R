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
#'
#' @importFrom BSgenome getSeq
#'
#' @return  fasta in specified directory.
#'
#'
#' @examples
#'
#' @export
#' @rdname outputTrees


outputTrees <- function(GDS, loci, dir, pops, nCores = 1, ploidy = 2, alleles = "seperate", minSites = 0.1, removeIndels = TRUE){


  store <- mclapply(1:length(loci), mc.cores = 4, function(locus){

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

      dna <- phyDat(genoMat)
      tree <- dist.ml(dna, model="JC69")
      tree <- upgma(tree)

      position <- seqGetData(gdsfile = GDS, var.name = "position")
      snpMid <- floor(median(position))
      nSites <- length(position)

      chr <- as.character(seqnames(locus))[[1]]
      start <- min(locus@ranges@start)
      end <- max(locus@ranges@start + width(locus))
      mid <- start+end / 2

      df <- data_frame(CHR = chr, Start = start, End = end, snpMid, nSites)

      if("Name" %in% colnames(locus@elementMetadata)){
        gNames <- paste(unique(locus$Name), collapse = ",")
        df <- bind_cols(df[1:3], data_frame(genes = gName), df[4:5])
      }

      df <- list(df, tree)

      df

    }


  })


  store <- unlist(store, recursive = FALSE)

  df <- lapply(seq(from = 1, to = length(store), by = 2), function(y){

    store[[y]]
  })

  df <- bind_rows(df)

  trees <- lapply(seq(from = 2, to = length(store), by = 2), function(y){

    store[[y]]
  })

  if(alleles == "seperate"){
    samples <- paste(rep(pops$Sample, each = ploidy), c(1:ploidy), sep = "/")
    pops <- rep(pops$Population, each = ploidy)
    pops <- data_frame(samples, pops)

  }

  class(trees) <- "multiPhylo"

  write.table(df, file = paste0(prefix, "_metaData.tsv"), quote = FALSE, row.names = FALSE)
  write.tree(trees, file = paste0(prefix, ".trees"))
  write.table(pops, file = paste0(prefix, "_groups.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)





}
