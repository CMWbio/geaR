% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getFeatures.R
\name{makeFeatures}
\alias{makeFeatures}
\title{getFeatures}
\usage{
makeFeatures(
  gffName,
  feature,
  nCores = 1,
  longestIsoform,
  includeRange,
  geneIDField,
  ...
)
}
\arguments{
\item{gffName}{\code{character} \cr
Path to gff to extract features from}

\item{feature}{\code{character}. \cr
Loci to import genotypes for. c("gene", "gene:exon", "gene:cds", "pseudogene", "lncRNA", "intergenic")}

\item{nCores}{\code{numeric} \cr
Number of cores to run in parallel}

\item{longestIsoform}{\code{logical} \cr
Only effects feature = "gene:exon". \cr
By default function will select the first entry for each gene. \cr
If \code{TRUE} is passed, will select the longest isoform for that gene. \cr
If there is a \code{"biotype"} field in the GFF it will only deal with those labeled as 'protein coding.' \cr
Selecting the longest isoform may be problematic if the \code{"biotype"} is missing.}

\item{includeRange}{\code{GRanges} or \code{GRangesList} \cr
Target ranges to extract features over.}

\item{contigMD}{\code{data.frame} \cr
Contig metadata. Contains two columns: \code{ID} with scaffold names, \code{length} with scaffold length. \cr
Used to select intergenic regions.}

\item{geneIdField}{\code{character} \cr
Name of the gene id field in the gff, \code{"Name"} by default.}
}
\value{
A \code{GRangesList} object where each element is a genomic feature
}
\description{
Creates feature loci from a gff file
}
\details{
Authour: Chris Ward \cr
Get genomic features such as 'genes', 'coding sequence', or 'exons' from gff.
}
\examples{

}
