#' @title The gear Object Class
#'
#' @description The gear Object Class
#'
#' @details This object contains \code{cog} objects to allow for multiple analyses to
#' be carried out on the same loci
#'
#' @return An object of class geaR

#' @rdname gear-class
#' @aliases gear-class
#' @export
setClass(
  "gear",
  slots = c(
    Loci = "GRangesList",
    Populations = "data.frame",
    Outgroup = "character",
    Args = "nullORargs",
    DiversityStatsFAST = "nullORdivFA",
    DiversityStatsFULL = "nullORdivFU",
    AdmixtureStats = "nullORadm",
    OutputLoci = "nullORoutL",
    OutputTrees = "nullORoutT"
  )
)


#### build empty class
.initializeGear <- function(){

  new("gear",
      Loci = GRangesList(GRanges("a:1-2")),
      Populations = data.frame(),
      Outgroup = character(),
      Args = .cogNull(),
      DiversityStatsFAST = .cogNull(),
      DiversityStatsFULL = .cogNull(),
      AdmixtureStats = .cogNull(),
      OutputLoci = .cogNull(),
      OutputTrees = .cogNull())



}

#' @title make gear object
#'
#' @description make a gear object from cogs
#'
#' @details This object contains \code{cog} objects to allow for multiple analyses to
#' be carried out on the same loci
#' 
#' @param loci \code{GRangesList} containing genomic loci to work on
#' @param populations \code{data.frame} contating two columns \code{Sample} and 
#' \code{Populaiton} 
#' @param outgroup \code{character} the population in the \code{populaitons} 
#' argument that corresponds to the outgroup
#' @param cogs \code{list} a list of cog objects for the analysis 
#'
#' @rdname makeGear
#' @export
makeGear <- function(loci, populations, outgroup = vector("character", 1), cogs = list()){

  gear <- .initializeGear()
  classes <- unlist(lapply(cogs, function(x) class(x)))
  stopifnot(length(unique(classes)) == length(classes))

  if(any(grepl(":", pops$Population))) stop("remove the ':' from population names ")
  
  gear@Loci <- loci
  gear@Populations <- populations
  gear@Outgroup <- outgroup

  if(any(classes == "cog.args")) gear@Args <- cogs[[which(classes == "cog.args")]]
  if(any(classes == "cog.diversityFAST")) gear@DiversityStatsFAST <- cogs[[which(classes == "cog.diversityFAST")]]
  if(any(classes == "cog.diversityFULL")) gear@DiversityStatsFULL <- cogs[[which(classes == "cog.diversityFULL")]]
  if(any(classes == "cog.admixture")) gear@AdmixtureStats <- cogs[[which(classes == "cog.admixture")]]
  if(any(classes == "cog.outputLoci")) gear@OutputLoci <- cogs[[which(classes == "cog.outputLoci")]]
  if(any(classes == "cog.outputTrees")) gear@OutputTrees <- cogs[[which(classes == "cog.outputTrees")]]

  return(gear)
}

