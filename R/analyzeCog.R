#' analyze cogs
#'
#' @description applies set methods to analyze a supplied cog
#'
#' @details 
#' backend methods to calculate cogs for a \code{gear} object
#'
#' @return the output of the supplied cog
#'
#'
#' @import furrr
#' @importFrom future plan
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#'
#'
#' @include cog-class.R
#' @rdname analyzeCog-methods
setGeneric("analyzeCog", function(cog, ...){
    standardGeneric("analyzeCog")
})

#' @aliases analyzeCog
setMethod("analyzeCog", signature(c(cog = "cog.NULL")),
          function(cog, ...){
              return(NULL)
          })

#' @aliases analyzeCog
setMethod("analyzeCog", signature(c(cog = "cog.diversityFULL")),
          function(cog, pairs, distMat,
                   arg, popList, seqname, start, end, windowMid, snpMid, nSites, locus, outgroup){
              
              
              #### build up stats to calculate
              if(cog@stats == "all") cog@stats <- c("pi", "dxy", "da", "dmin", "dmax", "Fst")
              
              
              #### set up dummy variables to build dataframe
              #### they will remain empty unless specified by cog@stats
              dxy <- c()
              pi <- c()
              da <- c()
              dmin <- c()
              dmax <- c()
              fst <- c()
              rndMin <- c()
              rndFeder <- c()
              rD <- c()
              
              #### conditionally calculates stats of interes based on cog@stats
              if("dxy" %in% cog@stats) dxy <- neisDxy(distMat, popList, pairs, ploidy = arg@ploidy)
              
              if("pi" %in% cog@stats) pi <- neisPi(distMat, popList, ploidy = arg@ploidy)
              
              if("dmin" %in% cog@stats) dmin <- dmin(distMat, popList, pairs, ploidy = arg@ploidy)
              
              if("dmax" %in% cog@stats) dmax <- dmax(distMat, popList, pairs, ploidy = arg@ploidy)
              
              if("da" %in% cog@stats) da <- neisDa(dxy, pi)
              
              if("Fst" %in% cog@stats) fst <- Nei82Fst(distMat, popList, pairs, ploidy = arg@ploidy, weighted = TRUE)
              
              if(!is.null(outgroup)){
                  if("RNDmin" %in% cog@stats) rndMin <- RND(istMat, popList, pairs, ploidy = arg@ploidy, outgroup = outgroup, type = "min")
                  
                  if("RNDfeder" %in% cog@stats) rndFeder <- RND(istMat, popList, pairs, ploidy = arg@ploidy, outgroup = outgroup, type = "feder")
                  
                  if("FTD" %in% cog@stats) rD <- relDxy(distMat, popList, ploidy = arg@ploidy, outgroup = outgroup)
                  
              }
              
              
              #### if a genes are used in construction of gear@Loci tehy will have a name
              #### retreive this and add this to the filename
              if("Name" %in% colnames(locus@elementMetadata)) {
                  gNames <- paste(unique(locus$Name), collapse = ",")
                  bind_cols(tibble(SeqName = seqname, Start = start, End = end, Gene = gNames, windowMid, snpMid, nSites), pi, dxy, da, dmin, dmax, fst, rndMin, rndFeder, rD)}
              else bind_cols(tibble(SeqName = seqname, Start = start, End = end, windowMid, snpMid, nSites), pi, dxy, da, dmin, dmax, fst, rndMin, rndFeder, rD)
          })


setMethod("analyzeCog", signature(c(cog = "cog.diversityFAST")),
          function(cog, pairs, distMat,
                   arg, popList, seqname, start, end, windowMid, snpMid, nSites, locus, outgroup, AF){
              
              
              
          })

#' @aliases analyzeCog
setMethod("analyzeCog", signature(c(cog = "cog.admixture")),
          function(cog, arg, pops,
                   locus, outgroup, AF){
              
              ### get the unique populations that do not belong to the outgroup 
              ### to construct comparisons from and to reduce the search space if 
              ### specified tests do not contain all populations in gear@Populaitons
              uTest <- unique(pops$Population)
              noOut <- uTest[uTest != outgroup]
              if(any(unlist(cog@fourPop) == "all")) {
                  combs <- combn(noOut, m = 3)
                  cog@fourPop <- lapply(seq(dim(combs)[[2]]), function(x){
                      x <- combs[,x]
                      c(x, outgroup)
                  })
              }
              
              p <- pops[pops$Population %in% unique(unlist(cog@fourPop)),]
              
              calcFourPop <- NULL
              if(length(AF)){
                  
                  calcFourPop <- map(cog@fourPop, function(x){
                      f4 <- .fourPop(AF, locus = locus, pops = pops, x = x)
                  })
                  scaf <- locus@seqnames[1]
                  st <- locus@ranges@start[1]
                  end <- st + sum(locus@ranges@width)
                  calcFourPop <- bind_cols(calcFourPop)
                  sM <- calcFourPop$snpMid
                  nS <- calcFourPop$nSites
                  
                  calcFourPop <- calcFourPop[!grepl("nSites", colnames(calcFourPop)) & !grepl("snpMid", colnames(calcFourPop))]
                  calcFourPop <- bind_cols(tibble(SeqName = as.character(scaf), Start = st, End = end, windowMid = (st + end)/2, snpMid = sM, nSites = nS), calcFourPop)
              } 
              
              return(calcFourPop)
              
          })

#' @aliases analyzeCog
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings writeXStringSet
setMethod("analyzeCog", signature(c(cog = "cog.outputLoci")),
          function(cog, genoMat, arg, locus, pops){
              
              #### set env variables
              ploidy <- arg@ploidy
              alleles <- cog@alleles
              removeIndels <- arg@removeIndels
              dir <- cog@outputDirectory
              
              #### fix the introduction of *s by gatk
              if(removeIndels){
                  genoMat[genoMat == "*"] <- NA
                  genoMat[complete.cases(genoMat),]
              }
              else genoMat[genoMat == "*"] <- "N"
              
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
                      pos <- rownames(genoMat)
                      
                      genoMat <- sapply(seq(from = 2, to = ncol(genoMat), by = 2), function(k){
                          
                          alleles <- paste(genoMat[,k-1], genoMat[,k], sep = "/")
                          alleles <- codevec[as.vector(alleles)]
                          
                      })
                      colnames(genoMat) <- unique(gsub("/.*", "", names))
                      rownames(genoMat) <- pos
                      
                  }
                  
                  if("gene" %in% colnames(locus@elementMetadata)){
                      
                      st <- min(locus@ranges@start)
                      ed <- min(end(locus))
                      
                      filename <- paste(as.character(GRanges(seqnames = locus@seqnames[1],
                                                             strand = locus@strand[1],
                                                             IRanges(start = st, end = ed))), collapse = ",")
                      
                      
                      filename <- paste0(locus$gene[[1]], "_", filename, ".fasta")
                      
                  }
                  
                  else{
                      
                      if("lociType" %in% colnames(locus@elementMetadata)){
                          filename <- paste0(as.character(locus@seqnames[1]), paste0(":", locus@elementMetadata$lociType[1]), round(mean(locus@ranges@start)))
                      }
                      
                      
                      filename <- paste0(filename, ".fasta")
                      
                  }
                  
                  genoMat <- t(genoMat)
                  genoMat <- split(genoMat, factor(rownames(genoMat), levels = rownames(genoMat)))
                  genoMat <- sapply(genoMat, paste, collapse="")
                  genoMat <- DNAStringSet(genoMat)
                  ref <- genoMat$Ref
                  genoMat <- genoMat[names(genoMat) != "Ref"]
                  
                  if(all(locus@strand == "-")) {
                      genoMat <- reverseComplement(genoMat)
                  }
                  genoMat <- c(DNAStringSet(ref), genoMat)
                  writeXStringSet(genoMat, paste0(dir, "/", filename))
                  
              }
              
              
          })


#' @aliases analyzeCog
#' @importFrom ape nj
#' @importFrom ape write.tree
setMethod("analyzeCog", signature(c(cog = "cog.outputTrees")),
          function(cog, GDS, arg, pops, locus, distMat){
              
              #### set env variables
              ploidy <- arg@ploidy
              alleles <- cog@alleles
              dir <- cog@outputDirectory
              
              
              
              
              if(length(distMat)){
                  
                  if(alleles == "concensus" & ploidy == 2){
                      
                      names <- colnames(distMat)
                      
                      distMat <- sapply(seq(from = 2, to = ncol(distMat), by = 2), function(k){
                          
                          sampleMean <- rowMeans(cbind(distMat[,k-1], distMat[,k]))
                          
                      })
                      
                      distMat <- sapply(seq(from = 2, to = nrow(distMat), by = 2), function(k){
                          
                          sampleMean <- rowMeans(cbind(distMat[k-1,], distMat[k,]))
                          
                      })
                      
                      
                      colnames(distMat) <- unique(gsub("/.*", "", names))
                      rownames(distMat) <- unique(gsub("/.*", "", names))
                      
                  }
                  
                  if(!any(is.na(distMat))){
                      tree <- nj(distMat)
                      
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
                      
                      
                      write.table(df, file = paste0(dir, "distanceTrees_metaData.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
                      write.tree(tree, file = paste0(dir, "distanceTrees.trees"), append = TRUE)
                      
                  }
              }
          })

