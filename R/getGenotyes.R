#' Import genotypes from GDS
#'
#' @description Imports genotypes across loci into matrix
#'
#' @details Authours: Chris Ward & Alastair Luddington
#' generates tiled (no step size) and sliding (step size) ranges for analysis with other functions.
#'
#'
#' @param GDS \code{GDS} object with variant data to import genotypes from
#' @param loci Loci to import genotypes for
#' @param minSites \code{numeric} minimum number of sites as a proportion of loci length
#' @param ncores \code{numeric} number of cores to run in parallel
#'
#' @return A \code{list} of \code{GRanges}
#'
#'
#' @import SeqArray
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges width
#' @importFrom car recode
#'
#' @examples
#'
#'
#'
#' @export
#' @rdname getGenotypes


getGenotypes <- function(GDS, locus, minSites = 0.5, nCores = 1, nuc = FALSE){

  stopifnot(minSites < 1 & minSites > 0)

  ## Iterating through the window list by scaffold
  by_scaffold <- lapply(names(loci), function(scaffoldName){

    ## Iterating through windows for current scaffold
    by_window <- pbmclapply(loci[[scaffoldName]], mc.cores = nCores, function(currentWindow){

      minSites <- minSites * width(currentWindow)

      ## Subsetting GDS file
      seqSetFilter(object = gds, currentWindow) ## Setting filter (i.e. window to get variants from)

      # read in genotypes and alleles
      genoArr <- seqGetData(gdsfile = GDS, var.name = "genotype")
      alleleArr <- seqGetData(gdsfile = GDS, var.name = "allele")

      ## Number of variants in window
      varNumber <- length(alleleArr)


      ## Filtering empty arrays/arrays with too few variants
      if(varNumber > minSites || is.null(arr)) {

        print("Incompatible window - not enough information")
        return(NULL)

      } else {


        ## Building window identifier
        chr <- as.vector(seqnames(currentWindow)@values) ## Getting the current scaffold
        # chr <- paste0("chr",chr)
        start <- start(currentWindow)
        end <- end(currentWindow)
        winID <- paste0(chr, ":", start, "..", end)

        ## Getting sample and variant IDs for col/row nameing
        samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
        # position <- seqGetData(gdsfile = GDS, var.name = "position")


        #Do I want to have a matrix of leave as an array?
        # get ploidy
        ploidy <- dim(genoArr)[1]
        ## Duplicating row-names (duce for splitting into haplotypes)
        #samples <- paste(rep(samples, each = ploidy), c(1:ploidy), sep = "/")


        if(nuc){
        # convert to nuceotides
        genoList <- lapply(1:varNumber, function(x){

          # get position, set sample names  and replace NA with N
          mat <- genoArr[,,x]
          mat[is.na(mat)] <- "N"

          alleles <- strsplit(alleleArr[x], ",")[[1]]
          geno <- 1:length(alleles) -1

          res <- matrix(alleles[match(mat, geno)], 2)
          mat <- ifelse(is.na(res), mat, res)
          mat <- c(mat)
          names(mat) <- paste(rep(samples, each = ploidy), c(1:ploidy), sep = "/")
          return(mat)
          })

        genoMat <- do.call(rbind, genoList)
        }
        else {

          genoMat <- t(apply(genoArr, MARGIN = 3, function(z){c(z)}))
          colnames(genoMat) <- paste(rep(samples, each = ploidy), c(1:ploidy), sep = "/")
          genoMat[is.na(genoMat)] <- "N"

        }

        geno <- union(genoMat, genoMat)

        dat <- lapply(geno, function(f){

          mat <- genoMat == f
          t(mat) %*% mat
        })

        dat <- Reduce("+", dat)

        distMat <- ( nrow(genoMat) - dat ) / nrow(genoMat)

      }

    })

  })

  ## Assigning scaffold names to list elements
  names(by_scaffold) <- names(windows_list)

  ## Returning named list of genotype matrices (split into haplotypes) by scaffold
  by_scaffold

}
