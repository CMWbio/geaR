#' Genotypes from GDS object
#'
#' @description Extracts genotype infromation from GDS objects and formats them for downstream analysis
#'
#' @details
#'
#'
#' @param fileName A \code{character} vector of length one containing the full path name for a Tabix indexed VCF
#'
#' @return A \code{tibble} containing populations statistics passed from \code{stat}
#'
#'
#' @import SeqArray
#'
#'
#' @examples
#'
#'
#' @export
#' @rdname GDSgenotypes


library(SeqArray)
## files for analysis
VCF <- "~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz"
GDS <- "~/Desktop/Tree-TipR/PlutellaSNP.GDS"

#get header from VCF will be needed for windows
header <- seqVCF_Header("~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz")$contig

# converty VCF to GDS
seqVCF2GDS(vcf.fn = "~/Desktop/Tree-TipR/PlutellaSNP.vcf.gz",
           parallel = 6, out.fn = "~/Desktop/Tree-TipR/PlutellaSNP.GDS")

## Parameters
windowSize <- 100000
step <- 100000
minSites <- 1000







#open GDS
gds <- seqOpen(GDS)




genotypes <-  seqGetData(gdsfile = gds, var.name = "genotype")

