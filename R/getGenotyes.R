#' Import genotypes from GDS
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
#' @rdname getGenotypes


getGenotypes <- function(GDS, locus, minSites = 0.5, nucleotide = FALSE){

  stopifnot(minSites < 1 & minSites > 0)

  minSites <- minSites * width(locus)

  ## Subsetting GDS file
  seqSetFilter(object = GDS, locus) ## Setting filter (i.e. window to get variants from)

  # read in genotypes and alleles
  genoArr <- seqGetData(gdsfile = GDS, var.name = "genotype")
  alleleArr <- seqGetData(gdsfile = GDS, var.name = "allele")

  ## Number of variants in window
  varNumber <- length(alleleArr)


  ## Filtering empty arrays/arrays with too few variants
  if(varNumber < minSites || is.null(genoArr)) {

    print("Incompatible window - not enough information")
    return(NULL)

  } else {


    ## Building window identifier
    chr <- as.vector(seqnames(locus)@values) ## Getting the current scaffold
    # chr <- paste0("chr",chr)
    start <- start(locus)
    end <- end(locus)
    winID <- paste0(chr, ":", start, "..", end)

    ## Getting sample and variant IDs for col/row nameing
    samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
    position <- seqGetData(gdsfile = GDS, var.name = "position")


    #Do I want to have a matrix of leave as an array?
    # get ploidy
    ploidy <- dim(genoArr)[1]
    ## Duplicating row-names (duce for splitting into haplotypes)
    #samples <- paste(rep(samples, each = ploidy), c(1:ploidy), sep = "/")


    if(nucleotide){
      # convert to nuceotides
      genoList <- lapply(1:varNumber, function(x){

        # get position and replace NA with N
        mat <- genoArr[,,x]
        mat[is.na(mat)] <- "N"

        # split allele string into vector
        alleles <- strsplit(alleleArr[x], ",")[[1]]
        # get genotype coding
        geno <- 1:length(alleles) -1

        # change RAW genotypes to nucleotide genotypes
        res <- matrix(alleles[match(mat, geno)], 2)
        mat <- ifelse(is.na(res), mat, res)
        # vectorize to order alleles and name
        mat <- c(mat)
        names(mat) <- paste(rep(samples, each = ploidy), c(1:ploidy), sep = "/")

        return(mat)
      })

      # bind all elements in list into matrix
      genoMat <- do.call(rbind, genoList)

    }
    else {

      # concatenate array across 3d (Variant) margin, set colnames then replace NA with N
      # replacing with N will convert all ints to chars
      genoMat <- t(apply(genoArr, MARGIN = 3, function(z){c(z)}))
      colnames(genoMat) <- paste(rep(samples, each = ploidy), c(1:ploidy), sep = "/")
      genoMat[is.na(genoMat)] <- "N"

    }

    #set variant positional information
    rownames(genoMat) <- position
    return(genoMat)
  }

}
