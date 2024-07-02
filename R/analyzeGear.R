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
#' @importFrom future multisession
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom magrittr extract2
#'
#' @rdname analyzeGear-methods
#' @export
setGeneric("analyzeGear",function(GDS, gear = NULL){
    standardGeneric("analyzeGear")
})

#' @aliases analyzeGear,character
#' @export
setMethod("analyzeGear", signature = "character",
          function(GDS, gear){
              GDS <- seqOpen(GDS, allow.duplicate = TRUE)
              analyzeGear(GDS, gear)
          })


#' @aliases analyzeGear,character
#' @export
setMethod("analyzeGear", signature = c(GDS = "SeqVarGDSClass"),
          function(GDS, gear){
              
              ### plan the future
              # plan(multisession, workers = gear@Args@nCores)
              
              ##### set up variables
              ploidy <- gear@Args@ploidy
              nCores <- gear@Args@nCores
              minSites <- gear@Args@minSites
              pairwiseDeletion <- gear@Args@pairwiseDeletion
              removeIndels <- gear@Args@removeIndels
              
              pops <- gear@Populations
              
              #### this is the bulk of the package, takes the loci object in 
              #### gear@Loci and runs each cog on it 
              data <- mclapply(seq_along(gear@Loci), mc.cores = nCores, function(x){
                  
                  ### set some dummy variables 
                  divFULL <- NULL
                  divFAST <- NULL
                  admix <- NULL
                  
                  
                  
                  locus <- gear@Loci[[x]]
                  
                  ### general locus metadata
                  seqname <- locus@seqnames@values[1]
                  start <- locus@ranges@start
                  end <- start + locus@ranges@width
                  start <- min(start)
                  end <- max(end)
                  windowMid <- (start + end) /2
                  
                  
                  #### test if it is necessary to read in actual genotypes
                  #### if the user just wants AF modules we dont need to read into mem
                  if(class(gear@DiversityStatsFULL) == "cog.diversityFULL" |
                     class(gear@OutputLoci) == "cog.outputLoci" | 
                     class(gear@OutputTrees) == "cog.outputTrees"){
                      
                      ### read in the raw genotype matrix
                      rawMat <- getGenotypes(GDS = GDS, locus = locus, minSites = minSites,
                                             raw = TRUE, ploidy = ploidy, pops = pops, removeIndels = removeIndels)
                      
                      if(length(rawMat)){
                          
                          ## read in genotype positional information
                          position <- seqGetData(gdsfile = GDS, var.name = "position")
                          
                          
                          ## snp metadata
                          snpMid <- floor(median(position))
                          nSites <- length(position)
                          if(length(position)){
                              
                              
                              ### get info from rawMat
                              ## genotype array
                              genoArr <- rawMat[[1]]
                              ## number of variants in array
                              varNumber <- rawMat[[2]]
                              ## sample names in array
                              samples <-rawMat[[3]]
                              
                              genoMat <- t(apply(genoArr, MARGIN = 3, function(z){c(z)}))
                              rownames(genoMat) <- position
                              colnames(genoMat) <- paste(rep(samples, each = gear@Args@ploidy), c(1:gear@Args@ploidy), sep = "/")
                              genoMat[is.na(genoMat)] <- "N"
                              
                              #### dummy variable to the distance matrix
                              distMat <- matrix()
                              
                              
                              #### only need to calc the distance matrix for "cog.diversityFULL" or
                              #### "cog.outputTrees"
                              needDist <- c(class(gear@DiversityStatsFULL),
                                            class(gear@OutputTrees)) %in% c("cog.diversityFULL",
                                                                            "cog.outputTrees")
                              
                              
                              ######### This chunk will resolve diversityFULL and outputTree modules
                              if(any(needDist)){
                                  distMat <- genoDist(genoMat, gear@Args@pairwiseDeletion)
                                  #### pregenerate populations list and pairs for pairwise dxy
                                  popList <- split(pops, pops$Population)
                                  pairs <- .generatePairs(popList)
                                  
                                  #### calculate diversity statistics, will return NULL if gear@DiversityStatsFULL
                                  #### slot is cog.NULL
                                  divFULL <- geaR:::analyzeCog(cog = gear@DiversityStatsFULL, pairs = pairs, distMat = distMat,
                                                        arg = gear@Args, popList, seqname, start, end, windowMid, snpMid, nSites, locus, outgroup = gear@Outgroup)
                                  #### output trees, will return NULL if gear@OutputTrees
                                  #### slot is cog.NULL
                                  analyzeCog(cog = gear@OutputTrees, GDS, arg = gear@Args, pops = gear@Populations, locus, distMat)
                              }
                              
                              ###### need to read in necleotides, distance matrix was slow as dick with freebayes MNPs
                              if(class(gear@OutputLoci) == "cog.outputLoci") {
                                  # convert genomat numeric genotypes 0123 to nuceotides ATCG
                                  genoMat <- .convertToNucleotide(GDS, varNumber, genoArr, removeIndels, ploidy, samples, position)
                                  
                                  #### output loci, will return NULL if gear@OutputTrees
                                  #### slot is cog.NULL
                                  analyzeCog(cog = gear@OutputLoci, genoMat, arg = gear@Args, locus,
                                             pops = gear@Populations, pos = position)
                              }
                          }
                          
                          
                          
                      }
                      
                  }
                  
                  ##### resolve diversityFAST and admixture modules
                  if(class(gear@DiversityStatsFAST) == "cog.diversityFAST" | 
                     class(gear@AdmixtureStats) == "cog.admixture") {
                      
                      # precompute allele frequencies
                      AF <- getAF(GDS, locus, pops = pops, minSites = minSites, refAllele = 0)
                      
                      # divFAST <- analyzeCog()
                      
                      #### calculate admixture statistics, will return NULL if gear@AdmixtureStats
                      #### slot is cog.NULL
                      admix <- geaR:::analyzeCog(cog = gear@AdmixtureStats, arg = gear@Args, pops = pops,
                                          locus, outgroup = gear@Outgroup, AF)
                      
                  }
                  
                  
                  ### this will save some memory if only one of diversity and admixture analyses are calculated
                  outList <- list(NULL, NULL)
                  if(!is.null(divFULL)) outList[[1]] <- divFULL 
                  if(!is.null(admix)) outList[[2]] <- admix
                  return(outList)  
                  
                  
                  
              })
              ## generate output
              if(class(gear@DiversityStatsFULL) == "cog.diversityFULL") gear@DiversityStatsFULL <- dplyr::bind_rows(lapply(data, extract2, 1))
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
