#' Calculates block jackkniffed mean for statistics
#'
#' @description summarize statistics using a block jacknife
#'
#' @details Authours: Chris Ward
#'
#'
#' @param x A \code{matrix} \cr
#' Allele genotypes for each individual
#' @param blockSize  \code{integer} default \code{1000000L}. The block size in
#' basepairs.
#' @param which \code{character,NULL} default \code{NULL}. Only calculates
#' across specified scaffolds/chromosomes. if  \code{NULL} will calculate
#' across all windows provided.
#' @param nCores \code{numeric} Number of cores allcated for calculations.
#' @param stats \code{character,NULL} default \code{NULL}. Only calculates
#' for specified statistics contained in the dataframe. if  \code{NULL} will
#' calculate across all statistics provided.
#' @param nCores \code{numeric} Number of cores allcated for calculations.
#' @return A \code{dataframe} with mean, jackSE, Zscore, median, and maximum for
#' each statistic.
#'
#' @rdname summarizeStats
#' @export


summarizeStats <- function(x, blockSize = 1000000L, which = NULL, nCores = 1, stats = NULL){

  cols <- colnames(x)

  if(!is.null(which)) x <- x[SeqName %in% which]


  scaffs <- as.character(unique(x$SeqName))

  xList <- split(x, as.character(x$SeqName))

  blocks <- lapply(1:length(xList), function(y){

    scaf <- xList[[y]]

    blockStart <- seq(min(scaf$Start), max(scaf$End), blockSize)

    blockEnd <- blockStart + blockSize


    tibble(chr = scaffs[y], blockStart, blockEnd)

  })

  blocks <- bind_rows(blocks)

  blocks <- split(blocks, 1:nrow(blocks))


  avgStats <- mclapply(seq_along(blocks), mc.cores = nCores, function(z){


    block <- blocks[[z]]

    avgStats  <-  x[
             x[["SeqName"]] == block[["chr"]] &
            x[["Start"]] >= block[["blockStart"]]&
            x[["End"]] <= block[["blockEnd"]],
            ]


    avgStats <- avgStats[,!cols %in% c("SeqName", "Start", "End", "windowMid", "snpMid", "Gene") &
                           ifelse(grepl(paste(paste0("_", stats), collapse = "|"), cols), TRUE, FALSE)]

    colMeans(avgStats, na.rm = TRUE)

  })

  avgStats <- do.call(rbind, avgStats)

  jackStats <- apply(avgStats, 2, function(k){

    u <- function(j) mean(j, na.rm = TRUE)

    jackVal <- bootstrap::jackknife(k, u)

    tibble(mean = mean(k, na.rm = TRUE), jackSE = jackVal$jack.se, Z = mean(k, na.rm = TRUE)/sd(k, na.rm = TRUE), median = median(k, na.rm = TRUE),
           minimum = min(k, na.rm = TRUE), maximum = max(k, na.rm = TRUE))

  })


  jackStatsDF <- bind_rows(jackStats)

  jackStatsDF["stat"] <- names(jackStats)

  jackStatsDF

}
