#' Make cogs 
#'
#' @description Builds the cog to pass to \code{makeGear}
#'
#' @details 
#' backend methods to build cogs for a \code{gear} object  \cr
#' Specify the \code{analysisType} and the corresponding argument set 
#'
#' @param analysisType \code{character}. Analysis type: \cr
#' @param ... See examples or the vignette for what arguments to pass in \code{...} \cr
#' \code{"args"}: Generic arguments for an analysis \cr
#' \code{"diversityFULL"}: Arguments spcific to the Diversity module \cr
#'
#' \itemize{
#' Currently supports: \cr
#' \item pi: Neis nucleotide diversity \cr
#' \item dxy: Neis absolute genetic distance \cr
#' \item da: Nei's ancesteral distance \cr
#' \item dmin: minimum hamming distance between samples \cr
#' \item dmax: maximum hamming distance between samples \cr
#' \item Fst: Nei (1982) \eqn{\gamma}st which estimates Fst \cr
#' \item RNDmin: minimum Relative node depth
#' \item federRND: Relative Node Depth using Feder (2005)}

#' \code{"admixture"}: Arguments spcific to the Admixture module \cr
#' \code{"outputLoci"}:  Arguments spcific to the Output Loci module \cr
#' \code{"outputTrees"}: Arguments spcific to the Output Trees module \cr
#'
#'
#' @return Output a \code{cog} object to 
#'
#' @examples
#' \dontrun{
#' 
#' # Argument cog
#' ## Necessary arguments are:
#' ## ploidy: The ploidy called in the genotyping process.
#' ## nCores: Number of cores to run analysis over.
#' ## minSites: The proportion of sites within the window with a genotype called. 
#' ## pairwiseDeletion: This will remove missing data in a pairwise manner i.e. will not count missing data as informative.
#' ## removeIndels: Should indels be removed from the analysis?
#' arg <- makeCog(analysisType = "args", ploidy = 2, nCores = 4, minSites = 0.2, pairwiseDeletion = TRUE, removeIndels = TRUE)
#' 
#' # Genetic Diversity cog
#' ## stats: What statistics should be run? 
#' diversity <- makeCog(analysisType = "diversityFULL", stats = "all")
#' 
#' # Admixture cog
#' ## threePop: structure for f3 stat 
#' ## fourPop: structure for f4 and fd stat 
#' admix <- makeCog(analysisType = "admixture", threePop = list(c("P1", "P2", "O")),
#'                  fourPop = list(c("P1", "P2", "P3", "O")))
#'                  
#' # Output Loci cog
#' ## outputDirectory: The directory to output fasta files to for downstream analysis
#' ## alleles: Output "seperate" or "collapsed" genotypes
#' outloci <-  makeCog(analysisType = "outputLoci", outputDirectory = "~/path/to/dir/", alleles = "seperate")
#' 
#' # Output Trees cog
#' ## outputDirectory: The directory to output fasta files to for downstream analysis
#' ## alleles: Output "seperate" or "collapsed" genotypes
#' outTrees <-  makeCog(analysisType = "outputTrees", outputDirectory = "~/path/to/dir/", alleles = "seperate")
#' }
#'
#' @rdname makeCogs-methods
#' @export
########### building the cogs classes
makeCog <- function(analysisType = "args", ...){

  l <- list(...)

  if(analysisType == "diversityFAST"){

    cog <- new("cog.diversityFAST", stats = l$stats)

    return(cog)
  }

  if(analysisType == "diversityFULL"){

    cog <- new("cog.diversityFULL", stats = l$stats)

    return(cog)
  }

  if(analysisType == "args"){

    cog <- new("cog.args", ploidy = l$ploidy, nCores = l$nCores, minSites = l$minSites, pairwiseDeletion = l$pairwiseDeletion, removeIndels = l$removeIndels)

    return(cog)
  }

  if(analysisType == "admixture"){

    cog <- new("cog.admixture", fourPop = l$fourPop, threePop = l$threePop)

    return(cog)
  }

  if(analysisType == "outputLoci"){

    cog <- new("cog.outputLoci", outputDirectory = l$outputDirectory, alleles = l$alleles, hapSamples = l$hapSamples)

    return(cog)
  }

  if(analysisType == "outputTrees"){

    cog <- new("cog.outputTrees", outputDirectory = l$outputDirectory, alleles = l$alleles)

    return(cog)
  }


}






