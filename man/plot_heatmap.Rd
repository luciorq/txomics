% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{plot_heatmap}
\alias{plot_heatmap}
\title{Plot a Heatmap From Expression Data}
\usage{
plot_heatmap(tx, sample_table, color_by = NULL, num = NULL,
  scale = "row", show_rownames = FALSE, ...)
}
\arguments{
\item{tx}{txi object, matrix or data frame with trancript abundance}

\item{sample_table}{metadata data frame, containing description
of each sample and experimental design}

\item{color_by}{variable from sample table to group samples in the plot}

\item{num}{number of genes to plot, filtering by variance between samples,
if num is negative, it is used the absolute num of genes with
smaller, if num is NULL all genes are used.
default = NULL.}

\item{scale}{logical, default = "row", normalize data by row, column or none}

\item{show_rownames}{logical, default = FALSE}

\item{...}{additional plot parameters, check pheatmap documentation for
details}
}
\description{
Plot a heatmap from expression data
}
