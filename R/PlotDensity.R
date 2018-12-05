#' Makes Density plots
#'
#' @description Plots stats into Density plots
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


ThreeTipDistance <- function(data, I1, I2, O){

  I1vO <- data[[paste0(O, "v", I1, "_dxy")]]
  I2vO <- data[[paste0(O, "v", I2, "_dxy")]]

  TTD <- tibble::data_frame(TTD = I2vO / (I1vO + I2vO))

  df <- dplyr::bind_cols(data[1:6], TTD)


  my_dens <- function(data, mapping, ...) {
    ggplot(data = data, mapping=mapping) +
      geom_density(..., alpha = 0.7, color = NA)
  }



  plot <- ggplot(df, aes(x = TTD)) +geom_density(aes(y = ..scaled..), colour = NA, fill = "red", alpha = 0.4) + xlim(0,1) + theme_minimal()

  ggplotly(plot)

  lowest <- filter(data, TTD < 0.3)

  TDout <- filter(df, SeqName %in% unique(lowest$SeqName))

 ggplot(TDout, aes(x = windowMid, y = TTD)) + geom_line() +
    facet_wrap(~SeqName, ncol = 1, scales = "free_x") + ylim (0, 0.7) + theme_linedraw()

  Outliers <- filter(data, SeqName %in% unique(lowest$SeqName))[c(1,4,10,12,14)]
  Outliers <- gather(Outliers, key = "Comp", value ="dxy", 3:5)

  Outliers <- dData
  Outliers <- Outliers[c(1,4,15,17,19,23,25)]
  Outliers <- gather(Outliers, key = "Comp", value ="dxy", 3:7)

  ggplot(Outliers, aes(x = windowMid, y = dxy, colour = Comp)) + geom_line() +
    facet_wrap(~SeqName, ncol = 1, scales = "free_x") + ylim (0, 0.7) + theme_linedraw()

  Outliers <- dData[c(1,4,7,8,9)]
  Outliers <- gather(Outliers, key = "Comp", value ="pi", 3:5)

  ggplot(Outliers, aes(x = windowMid, y = pi, colour = Comp)) + geom_line() +
    facet_wrap(~SeqName, ncol = 1, scales = "free_x") + theme_linedraw()




}
