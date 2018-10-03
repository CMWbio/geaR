#' Combines gffs output by BUSCO
#'
#' @description Combines gffs from BUSCO annotation of a genome
#'
#' @details Authour: Chris Ward
#'
#' @param inDir directory to input single gffs from. gffs should be in *geneName*.gff nomenclature
#' @param outDir directory to output combined gff
#' @param geneList only import genes from this vector/list
#'
#' @importFrom rtracklayer import.gff
#'
#' @return A gff contatining all genes in geneList
#'
#'
#' @examples
#'
#'
#'
#' @export
#' @rdname combineBUSCOgff

combineBUSCOgff <- function(inDir, outDir, outName, geneList = NULL){

  if(length(geneList)){

    names <- paste0(inDir, "/", geneList, ".gff")

  } else {
    names <- list.files(inDir, "*.gff", full.names = TRUE)
  }

  grList <- lapply(names, function(x){

    gr <- import.gff(x)
    gene <- gsub(inDir, "", x)
    gene <- gene <- gsub("/", "", gene)
    gene <- gene <- gsub("\\.gff", "", gene)
    gr <- as_data_frame(gr)
    gr["ID"] <- gene
    gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = TRUE)

  })

  gr <- do.call("c", grList)
  export.gff3(gr, paste0(outDir, "/", outName, "gff"))
}
