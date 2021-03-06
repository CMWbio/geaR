% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGenotypes.R
\name{getGenotypes}
\alias{getGenotypes}
\title{Import genotypes from GDS}
\usage{
getGenotypes(
  GDS,
  locus = NULL,
  minSites = 0.5,
  raw = FALSE,
  ploidy = 2,
  pops = NULL,
  removeIndels = TRUE
)
}
\arguments{
\item{GDS}{\code{GDS} object with variant data to import genotypes from}

\item{locus}{\code{GRanges} Locus to import genotypes for}

\item{minSites}{\code{numeric} minimum number of sites as a proportion of loci length}

\item{ploidy}{\code{numeric} ploidy of sample}

\item{pops}{\code{tibble} populaiton dataFrame}

\item{removeIndels}{removes indels}

\item{nucleotide}{\code{logical} Import RAW genotypes or nucleotides}
}
\value{
A \code{matrix} of genotypes
}
\description{
Imports genotypes across a locus into matrix
}
\details{
Authours: Chris Ward & Alastair Ludington
Uses a GRanges locus to import genotypes, either nucleotide or RAW, from a GDS file
}
\examples{



}
