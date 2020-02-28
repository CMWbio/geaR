#' @importFrom methods as

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

    cog <- new("cog.args", ploidy = l$ploidy, nCores = l$nCores, minSites = l$minSites, pairwiseDeletion = l$pairwiseDeletion)

    return(cog)
  }

  if(analysisType == "admixture"){

    cog <- new("cog.admixture", fourPop = l$fourPop, threePop = l$threePop)

    return(cog)
  }

  if(analysisType == "outputLoci"){

    cog <- new("cog.outputLoci", outputDirectory = l$outputDirectory, alleles = l$alleles, removeIndels = l$removeIndels)

    return(cog)
  }

  if(analysisType == "outputTrees"){

    cog <- new("cog.outputTrees", outputDirectory = l$outputDirectory, alleles = l$alleles, removeIndels = l$removeIndels)

    return(cog)
  }


}








