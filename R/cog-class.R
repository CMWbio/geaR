#' @title Cog classes
#'
#' @description Cog classes
#'
#' @details cog classes to fit into gear object slots
#'
#' @return An object of class cog.X
#'

#' @rdname cog-class
#' @aliases cog-class,cog.args
#' @export
setClass(
  "cog.args",
  slots = c(
    ploidy = "numeric",
    nCores = "numeric",
    minSites = "numeric",
    pairwiseDeletion = "logical",
    removeIndels = "logical"
  )
)

#' @rdname cog-class
#' @aliases cog-class,cog.diversityFAST
#' @export
setClass(
  "cog.diversityFAST",
  slots = c(
    stats = "character"
  )
)

#' @rdname cog-class
#' @aliases cog-class,cog.diversityFULL
#' @export
setClass(
  "cog.diversityFULL",
  slots = c(
    stats = "character"
  )
)

#' @rdname cog-class
#' @aliases cog-class,cog.admixture
#' @export
setClass(
  "cog.admixture",
  slots = c(
    fourPop = "list",
    threePop = "list"
    )
)

#' @rdname cog-class
#' @aliases cog-class,cog.outputLoci
#' @export
setClass(
  "cog.outputLoci",
  slots = c(
    outputDirectory = "character",
    alleles = "character",
    removeIndels = "logical"
  )
)


#' @rdname cog-class
#' @aliases cog-class,cog.outputTrees
setClass(
  "cog.outputTrees",
  slots = c(
    outputDirectory = "character",
    alleles = "character",
    removeIndels = "logical"
  )
)


#' @rdname cog-class
#' @aliases cog-class,cog.NULL
#' @export
setClass(
  "cog.NULL",
  slots = c(
    null = "NULL"
  )
)



#### build empty cogs

.cogNull <- function() new("cog.NULL")



#### set unions


setClassUnion("nullORargs", c("cog.args", "cog.NULL"))
setClassUnion("nullORdivFU", c("cog.diversityFULL", "cog.NULL", "data.frame"))
setClassUnion("nullORdivFA", c("cog.diversityFAST", "cog.NULL", "data.frame"))
setClassUnion("nullORadm", c("cog.admixture", "cog.NULL", "data.frame"))
setClassUnion("nullORoutL", c("cog.outputLoci", "cog.NULL", "character"))
setClassUnion("nullORoutT", c("cog.outputTrees", "cog.NULL", "character"))




##### building empty

.cogArgs <- function() makeCog(analysisType = "args", ploidy = numeric(), nCores = numeric(), minSites = numeric(), pairwiseDeletion = logical(), removeIndels = logical())

.cogDivFAST <- function() makeCog(analysisType = "diversityFAST", stats = character())

.cogDivFULL <- function() makeCog(analysisType = "diversityFULL", stats = character())

.cogAd <- function() makeCog(analysisType = "admixture", fourPop = list(), threePop = list())

.cogOutLoc <- function() makeCog(analysisType = "outputLoci", outputDirectory = character(), alleles = character(),
                                 removeIndels = logical())

.cogOutTre <- function() makeCog(analysisType = "outputTrees", outputDirectory = character(), alleles = character())


