#' analyzeGear
#'
#' @description analyzes the object class \code{gear}
#'
#' @details 
#' analyszes a supplied GDS using a \code{gear} object containing the desired 
#' cogs.
#' 
#' @param GDS \code{character} contating the path to a GDS file 
#' or an already imported \code{SeqVarGDSClass} 
#' @param gear \code{gear} an object created with \code{makeGear()} contating the desired analyses. \cr
#' Analyses that can be carried out are (see \code{?makeCog() for details}): \cr
#' Genetic Diversity
#' Admixture
#' Output Loci 
#' Output Trees
#'
#' @return the supplied object with cogs replaced by results
#'
#' @import furrr
#' @importFrom future plan
#' @importFrom future multiprocess
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom magrittr extract2
#'
#' @rdname analyzeGear-methods
#' @export
setGeneric("analyzeGear",function(GDS, gear = NULL, ...){
    standardGeneric("analyzeGear")
})

#' @aliases getDiversityStats,character
#' @export
setMethod("analyzeGear", signature = "character",
          function(GDS, ...){
              GDS <- seqOpen(GDS, allow.duplicate = TRUE)
              getDiversityStats(GDS, ...)
          })


#' @aliases getDiversityStats,character
#' @export
setMethod("analyzeGear", signature = c(GDS = "SeqVarGDSClass"),
          function(GDS, ...){
              
              ### plan the future
              plan(multiprocess, workers = gear@Args@nCores)
              
              #### this is the bulk of the package, takes the loci object in 
              #### gear@Loci and runs each cog on it 
              data <- .quiet(future_map(seq_along(gear@Loci), function(x, gear, Gname){
                  
                  ##### set up variables
                  ploidy <- gear@Args@ploidy
                  nCores <- gear@Args@nCores
                  minSites <- gear@Args@minSites
                  pairwiseDeletion <- gear@Args@pairwiseDeletion
                  removeIndels <- gear@Args@removeIndels
                  
                  pops <- gear@Populations
                  
                  GDS2 <- seqOpen(gds.fn = Gname, allow.duplicate = TRUE)
                  
                  locus <- gear@Loci[[x]]
                  
                  ### read in the raw genotype matrix to convert to nucleotide if 
                  ### gear@OutputLoci is an analysis
                  
                  rawMat <- getGenotypes(GDS = GDS2, locus = locus, minSites = minSites,
                                         raw = TRUE, ploidy = ploidy, pops = pops, removeIndels = removeIndels)
                  
                  ### get info from rawMat
                  ## genotype array
                  genoArr <- rawMat[[1]]
                  ## number of variants in array
                  varNumber <- rawMat[[2]]
                  ## sample names in array
                  samples <-rawMat[[3]]
                  
                  ### get indexed genotype positions
                  position <- seqGetData(gdsfile = GDS2, var.name = "position")
                  
                  ### convert from raw geno array to genotype matrix in order to process downstream
                  genoMat <- t(apply(genoArr, MARGIN = 3, function(z){c(z)}))
                  rownames(genoMat) <- position
                  colnames(genoMat) <- paste(rep(samples, each = gear@Args@ploidy), c(1:gear@Args@ploidy), sep = "/")
                  genoMat[is.na(genoMat)] <- "N"
                  
                  ### general locus metadata
                  seqname <- locus@seqnames@values[1]
                  start <- locus@ranges@start
                  end <- start + locus@ranges@width
                  start <- min(start)
                  end <- max(end)
                  windowMid <- (start + end) /2
                  snpMid <- floor(median(position))
                  nSites <- length(position)
                  
                  if(length(genoMat)){
                      
                      #### dummy variable to the distance matrix
                      distMat <- matrix()
                      
                      
                      #### only need to calc the distance matrix for "cog.diversityFULL" or
                      #### "cog.outputTrees"
                      needDist <- c(class(gear@DiversityStatsFULL),
                                    class(gear@OutputTrees)) %in% c("cog.diversityFULL",
                                                                    "cog.outputTrees")
                      
                      if(any(needDist)){
                          distMat <- genoDist(genoMat, gear@Args@pairwiseDeletion)
                          #### pregenerate populations list and pairs for pairwise dxy
                          popList <- split(pops, pops$Population)
                          pairs <- .generatePairs(popList)
                          
                          #### calculate diversity statistics, will return NULL if gear@DiversityStatsFULL
                          #### slot is cog.NULL
                          divFULL <- analyzeCog(cog = gear@DiversityStatsFULL, pairs = pairs, distMat = distMat,
                                                arg = gear@Args, popList, seqname, start, end, windowMid, snpMid, nSites, locus, outgroup = gear@Outgroup)
                          #### output trees, will return NULL if gear@OutputTrees
                          #### slot is cog.NULL
                          analyzeCog(cog = gear@OutputTrees, GDS2, arg = gear@Args, pops = gear@Populations, locus, distMat)
                      }
                      
                      #divFAST <- analyzeCog(cog = gear@DiversityStatsFAST)
                      
                      #### calculate admixture statistics, will return NULL if gear@AdmixtureStats
                      #### slot is cog.NULL
                      admix <- analyzeCog(cog = gear@AdmixtureStats, arg = gear@Args, pops = pops,
                                          locus, outgroup = gear@Outgroup, GDS2)
                      
                      ###### need to read in necleotides, distance matrix was slow as dick with freebayes MNPs
                      if(class(gear@OutputLoci) == "cog.outputLoci") {
                          # convert genomat numeric genotypes 0123 to nuceotides ATCG
                          genoMat <- .convertToNucleotide(GDS2, varNumber, genoArr, removeIndels, ploidy, samples, position)
                          
                          #### output loci, will return NULL if gear@OutputTrees
                          #### slot is cog.NULL
                          analyzeCog(cog = gear@OutputLoci, genoMat, arg = gear@Args, locus,
                                     pops = gear@Populations)
                      }
                  }
                  ### this will save some memory if only one of diversity and admixture analyses are calculated
                  outList <- vector(mode = "list", length = sum(c(!is.null(divFULL), !is.null(admix))))
                  if(!is.null(divFULL)) outList[[1]] <- divFULL 
                  if(!is.null(admix)) outList[[2]] <- admix
                  seqClose(GDS2)
                  return(outList)
              }, gear = gear, Gname = GDS$filename))
              plan(sequential)
              ## generate output
              if(class(gear@DiversityStatsFULL) == "cog.diversityFULL") gear@DiversityStatsFULL <- bind_rows(lapply(data, extract2, 1))
              if(class(gear@AdmixtureStats) == "cog.admixture") gear@AdmixtureStats <- bind_rows(lapply(data, extract2, 2))
              if(class(gear@OutputLoci) == "cog.outLoci") gear@OutputLoci <- paste0("Loci output to \"", gear@OutputLoci@outputDirectory)
              if(class(gear@OutputTrees) == "cog.outTrees") gear@OutputTrees <- paste0("Trees output to \"", gear@OutputTrees@outputDirectory)
              
              return(gear)
          })


.generatePairs <- function(popList, ...){
    nam <- names(popList)
    pairs <- outer(nam, nam, paste, sep = "///")
    pairs <- lapply(pairs,function(z){
        pair <- unlist(strsplit(z, split = "///"))
        if(pair[1] == pair[2]) pair <- NULL
        pair <- sort(pair)
        return(pair)
    })
    pairs <- pairs[!duplicated(pairs)]
    Filter(Negate(is.null), pairs)
}


.convertToNucleotide <- function(GDS, varNumber, genoArr, removeIndels, ploidy, samples, position){
    
    alleleArr <- seqGetData(gdsfile = GDS, var.name = "allele")
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
        
        ## remove alleles of unequal lengths ie insertions or deleletions
        if(removeIndels & !length(unique(nchar(mat))) == 1) return(NULL)
        
        
        mat <- matrix(mat, nrow = 1)
        rownames(mat) <- position[x]
        colnames(mat) <- paste(rep(samples, each = ploidy), c(1:ploidy), sep = "/")
        mat
    })
    
    # bind all elements in list into matrix
    do.call(rbind, genoList)
    
}

### hadley http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
.quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
} 