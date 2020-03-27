#' getDivFromGear
#'
#' @description Extract Diversity metrics from a gear object
#' 
#' @param gear \code{gear} object user has analyzed
#'
#' @return dataframe of results
#'
#' @rdname getDivFromGear
#' @export

getDivFromGear <- function(gear){
    return(gear@DiversityStatsFULL)
}

#' getAdmixFromGear
#'
#' @description Extract admixture metrics from a gear object
#' 
#' @param gear \code{gear} object user has analyzed
#'
#' @return dataframe of results
#'
#' @rdname getAdmixFromGear
#' @export

getAdmixFromGear <- function(gear){
    return(gear@AdmixtureStats)
}
