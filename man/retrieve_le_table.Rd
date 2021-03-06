% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{retrieve_le_table}
\alias{retrieve_le_table}
\title{Leading Edge Genes By Gene Set}
\usage{
retrieve_le_table(..., names = NULL)
}
\arguments{
\item{...}{gsea_analysis results data frames}

\item{names}{vector with names to be used for each group}
}
\description{
Retrieve table with leading edge genes by gene set, for multiple
experiments the column occurrence will have the number of hits
}
