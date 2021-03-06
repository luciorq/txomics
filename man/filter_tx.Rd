% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{filter_tx}
\alias{filter_tx}
\title{Filter Imported Transcripts}
\usage{
filter_tx(tx, sample_table, var_column = NULL, var_levels = NULL)
}
\arguments{
\item{tx}{txi object}

\item{sample_table}{metadata data frame, containing description
of each sample and experimental design}

\item{var_column}{variable to filter the input}

\item{var_levels}{levels of var_column to filter the input}
}
\value{
A \code{txi} object
}
\description{
Filter salmon quantification files
in a \code{txi} format
based on metadata column and levels
}
