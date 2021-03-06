% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{le_analysis}
\alias{le_analysis}
\title{Perform Leading Edge Analysis}
\usage{
le_analysis(gsea_res, set)
}
\arguments{
\item{gsea_res}{result from gsea_analysis}

\item{set}{gene set to be extracted}
}
\value{
A \code{vector} of leading edge genes for the set
}
\description{
Extract Leading Edge genes from GSEA results.
}
