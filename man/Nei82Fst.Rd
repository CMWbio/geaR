% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weightedFst.R
\name{Nei82Fst}
\alias{Nei82Fst}
\title{Fst from hamming distance}
\usage{
Nei82Fst(distMat, popList, pairs, ploidy, weighted = TRUE)
}
\arguments{
\item{distMat}{hamming distance matrix calculated using the hidden function \code{genoDist}}

\item{popList}{List of populations made from \code{pops} dataframe provided by users}

\item{pairs}{list of populaiton pairs generated in \code{getDiversityStats}}

\item{ploidy}{\code{numeric} number of chromosomes}

\item{weighted}{\code{logical} \cr
Weight by sample subpopulaiton size.
If \code{TRUE} mean nuceotide diversity within the subpopulaiton will be weigthed by sample size. \cr
With weighting 6 samples in p1 and 4 in p2, will result in 0.6*pi_p1 + 0.4*pi_p2 \cr
Alternatively, without weighting each subpopulaiton is assumed to contribute equally to the mean, (pi_p1 + pi_p2)/2}
}
\value{
A \code{dataframe} Fst between two populations
}
\description{
Calculates sample size weighted or unweighted Fst-like \eqn{\gamma}st from Nei (1982)
}
\details{
Authours: Chris Ward Uses the mean nucleotide diversity within each subpopulaiton to calculate subpopulaiton pi and divides by the total nucleotide diversity.
}
