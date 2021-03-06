% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTopSigGenes.R
\docType{methods}
\name{getTopSigGenes}
\alias{getTopSigGenes}
\alias{getTopSigGenes,BioData-method}
\title{obtain a list of most significant genes per comparison froma Seurat or Cpp stats test.}
\usage{
getTopSigGenes(x, name, n = 100)
}
\arguments{
\item{x}{the BioData object}

\item{name}{the grouping name (samples column name)}

\item{n}{the amount of genes per stat test patrt (default 100)}
}
\description{
obtain a list of most significant genes per comparison froma Seurat or Cpp stats test.
}
