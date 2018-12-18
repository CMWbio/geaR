#' get codon positions from fasta
#'
#' @description Imports genotypes across a locus into matrix
#'
#' @details Authours: Chris Ward & Alastair Ludington
#' Uses a GRanges locus to import genotypes, either nucleotide or RAW, from a GDS file
#'
#'
#' @param GDS \code{GDS} object with variant data to import genotypes from
#' @param locus \code{GRanges} Locus to import genotypes for
#' @param minSites \code{numeric} minimum number of sites as a proportion of loci length
#' @param nucleotide \code{logical} Import RAW genotypes or nucleotides
#'
#'
#' @return A \code{matrix} of genotypes
#'
#'
#' @import SeqArray
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges width
#' @examples
#'
#'
#'
#' @export
#' @rdname getCodonsFromFasta



getCodonsFromFasta <- function(genome = genome, exons){

    genes <- getSeq(genome, exons)

  Catgenes <- lapply(1:length(genes), function(x){

    unlist(genes[[x]])


   })

  names(Catgenes) <- names(genes)


ct <- c()
 codList <- lapply(1:length(Catgenes), function(x){

   ct <<-x

  if((length(Catgenes[[x]])/3)%%1 == 0 && !grepl("N", as.character(Catgenes[[x]]), ignore.case = TRUE)){

  dnaCodons <- codons(Catgenes[[x]]) %>% as.character() %>% tolower()

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
              met = "aug",
              Phe = "tt[tc]",
              Trp = "tgg",
              Tyr = "ta[tc]",
              Stp = "ta[ag]|tga"
  )

  FourFoldlookUP <- lookUP[1:8]


  match <- t(sapply(FourFoldlookUP, grepl, dnaCodons))

  redMatch <- Reduce(f = "|", split(match, rownames(match)))

  dnaCodons4fold <-  dnaCodons[redMatch]
  count <- table(dnaCodons4fold)

  RSCU <- lapply(FourFoldlookUP, function(z){

    countSplit <- count[grepl(z, rownames(count))]

    if(sum(countSplit) != 1) {
      RSCU <- lapply(countSplit, function(c){

        (4*c)/sum(countSplit)

      }) %>% unlist()

      #any(RSCU < 0.9 | RSCU > 1.4)
    }
    else NULL
  }) %>% unlist()

  df <- data_frame(codon = names(RSCU), RSCU = RSCU) %>% separate(codon, c("AA", "codon"))

  aaCodons <- dnaCodons

  s <- lapply(1:length(lookUP), function(z){

    aaCodons[grepl(lookUP[z], aaCodons)] <<- names(lookUP)[z]

  })

  rm(s)
  frst <- 1:length(dnaCodons)

  codonDF <- data_frame(AA = aaCodons, codon = dnaCodons, frst = 1:length(dnaCodons)*3-2, scnd = 1:length(dnaCodons)*3-1, thrd = 1:length(dnaCodons)*3)

  codonList <- list(RCSU = df, codonDF = codonDF)



  }


  })

 names(codList) <- names(FL)
 codList <- Filter(Negate(is.null), codList)


}
