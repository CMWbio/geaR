#' Calculate diversity statistics
#'
#' @description Imports genotypes and calculates diversity formula
#'
#' @details Authour: Chris Ward
#' Uses a \code{list} containing \code{GRanges} or a \code{GRangesList} to import genotypes. Diversity statistics are then calculated.
#'
#' @param GDS \code{GDS} object with variant data to import genotypes from.
#' @param loci \code{GRangesList}. Loci to import genotypes for.
#' @param minSites \code{numeric} minimum number of sites as a proportion of loci length. Default 0.5 ie 50 percent
#' @param nCores \code{numeric} number of cores to run in parallel
#' @param pops \code{data_frame} containing Sample ID and poplations. Default is to treat individuals as populations \code{c("none")}
#' @param stats \code{character}. Vector containing all diversity stats to calculate. default \code{c("all")}
#' @param ploidy \code{numeric} number of alleles
#' @param pairwiseDeletion \code{numeric} nshould Ns be removed from pairwise comparisons
#'
#'
#' @return A \code{data_frame} of selected Diversity statistics
#'
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#'
#' @rdname getDiversityStats-methods
#' @export



setGeneric("getDiversityStats",function(GDS, loci, minSites = 0.5, nCores = 1, pops = NULL, stats = "all", ploidy = 2, pairwiseDeletion, ...){
  standardGeneric("getDiversityStats")
})

#' @aliases getDiversityStats,character
#' @rdname getDiversityStats-methods
#' @export
setMethod("getDiversityStats", signature = "character",
          function(GDS, ...){
            GDS <- seqOpen(GDS, allow.duplicate = TRUE)
            getDiversityStats(GDS, ...)
          })


#' @aliases getDiversityStats,SeqVarsGDSClass-GRangesList
#' @rdname getDiversityStats-methods
#' @export
setMethod("getDiversityStats", signature(GDS = "SeqVarGDSClass",
                                         loci = "GRangesList"),
          function(GDS, loci, minSites = 0.5, nCores = 1, pops = NULL, stats = "all", ploidy = 2, pairwiseDeletion){


            if(stats == "all") stats <- c("pi", "dxy", "da", "dmin", "dmax", "Fst")

            if(length(pops) == 0){
              samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
              pops <- data_frame(Sample = samples, Population = samples)
            }

            popList <- split(pops, pops$Population)

            nam <- names(popList)
            pairs <- outer(nam, nam, paste, sep = "///")
            pairs <- lapply(pairs,function(z){
              pair <- unlist(strsplit(z, split = "///"))
              if(pair[1] == pair[2]) pair <- NULL
              pair <- sort(pair)
              return(pair)
            })
            pairs <- pairs[!duplicated(pairs)]
            pairs <- Filter(Negate(is.null), pairs)


            stopifnot(minSites < 1 & minSites > 0)
            ## Iterating through the window list by scaffold
            div <- pbmclapply(seq(length(loci)), mc.cores = nCores, function(locusN){
              locus <- loci[[locusN]]

              genoMat <- getGenotypes(GDS = GDS, locus = locus, minSites = minSites, nucleotide = FALSE, ploidy = ploidy, pops = pops)

              seqname <- locus@seqnames@values[1]
              start <- locus@ranges@start
              end <- start + locus@ranges@width

              start <- min(start)
              end <- max(end)
              windowMid <- (start + end) /2



              dxy <- c()
              pi <- c()
              da <- c()
              dmin <- c()
              dmax <- c()
              fst <- c()


              position <- seqGetData(gdsfile = GDS, var.name = "position")
              snpMid <- floor(median(position))
              nSites <- length(position)

              if(length(genoMat)){

                distMat <- geaR::genoDist(genoMat, pairwiseDeletion)

                if("dxy" %in% stats) dxy <- neisDxy(distMat, popList, pairs, ploidy = ploidy)

                if("pi" %in% stats) pi <- neisPi(distMat, popList, ploidy = ploidy)

                if("dmin" %in% stats) dmin <- dmin(distMat, popList, pairs, ploidy = ploidy)

                if("dmax" %in% stats) dmax <- dmax(distMat, popList, pairs, ploidy = ploidy)

                if("da" %in% stats) da <- neisDa(dxy, pi)

                if("fst" %in% stats) weightedFst(distMat, popList, pairs, ploidy = ploidy)

              }




              if("gene" %in% colnames(locus@elementMetadata)) {
                gNames <- paste(unique(locus$gene), collapse = ",")
                bind_cols(data_frame(SeqName = seqname, Start = start, End = end, Gene = gNames, windowMid, snpMid, nSites), pi, dxy, da, dmin, dmax, fst)}
              else bind_cols(data_frame(SeqName = seqname, Start = start, End = end, windowMid, snpMid, nSites), pi, dxy, da, dmin, dmax, fst)



            })

            dplyr::bind_rows(div)

          })


#' @aliases getDiversityStats,SeqVarsGDSClass-GRanges
#' @rdname getDiversityStats-methods
#' @export
setMethod("getDiversityStats", signature(GDS = "SeqVarGDSClass",
                                         loci = "GRanges"),
          function(GDS, loci, minSites = 0.5, nCores = 1, pops = NULL, stats = "all", ploidy = 2, pairwiseDeletion){

            if(stats == "all") stats <- c("pi", "dxy", "da", "dmin", "dmax", "Fst")

            if(length(pops) == 0){
              samples <- seqGetData(gdsfile = GDS, var.name = "sample.id")
              pops <- data_frame(Sample = samples, Population = samples)
            }

            popList <- split(pops, pops$Population)

            nam <- names(popList)
            pairs <- outer(nam, nam, paste, sep = "///")
            pairs <- lapply(pairs,function(z){
              pair <- unlist(strsplit(z, split = "///"))
              if(pair[1] == pair[2]) pair <- NULL
              pair <- sort(pair)
              return(pair)
            })
            pairs <- pairs[!duplicated(pairs)]
            pairs <- Filter(Negate(is.null), pairs)


            stopifnot(minSites < 1 & minSites > 0)
            ## Iterating through the window list by scaffold

            genoMat <- getGenotypes(GDS = GDS, locus = loci, minSites = minSites, nucleotide = FALSE, ploidy = ploidy, pops = pops)

            seqname <- loci@seqnames@values[1]
            start <- loci@ranges@start
            end <- start + loci@ranges@width

            start <- min(start)
            end <- max(end)
            windowMid <- (start + end) /2


            dxy <- c()
            pi <- c()
            da <- c()
            dmin <- c()
            dmax <- c()
            fst <- c()


            if(length(genoMat)){

              distMat <- geaR::genoDist(genoMat, pairwiseDeletion)

              if("dxy" %in% stats) dxy <- neisDxy(distMat, popList, pairs, ploidy = ploidy)

              if("pi" %in% stats) pi <- neisPi(distMat, popList, ploidy = ploidy)

              if("dmin" %in% stats) dmin <- dmin(distMat, popList, pairs, ploidy = ploidy)

              if("dmax" %in% stats) dmax <- dmax(distMat, popList, pairs, ploidy = ploidy)

              if("da" %in% stats) da <- neisDa(dxy, pi)

              if("fst" %in% stats) weightedFst(distMat, popList, pairs, ploidy = ploidy)

              snpMid <- floor(mean(as.numeric(rownames(genoMat)), na.rm = TRUE))
              nSites <- nrow(genoMat)

            }
            else {
              position <- seqGetData(gdsfile = GDS, var.name = "position")
              snpMid <- floor(median(position))
              nSites <- length(position)
            }




            if("gene" %in% colnames(locus@elementMetadata)) {
              gNames <- paste(unique(locus$gene), collapse = ",")
              bind_cols(data_frame(SeqName = seqname, Start = start, End = end, Gene = gNames, windowMid, snpMid, nSites), pi, dxy, da, dmin, dmax, fst)}
            else bind_cols(data_frame(SeqName = seqname, Start = start, End = end, windowMid, snpMid, nSites), pi, dxy, da, dmin, dmax, fst)



          })

