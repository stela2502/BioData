% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_Seurat.R
\docType{methods}
\name{as_Seurat}
\alias{as_Seurat}
\alias{as_Seurat,BioData-method}
\title{convert BioData to Seurat}
\usage{
as_Seurat(x, group = NULL, fromRaw = T, norm = TRUE)
}
\arguments{
\item{x}{The BioData object}

\item{group}{Which group should be used as identity in the seurat object (default=NULL)}

\item{fromRaw}{create the Seurat object from the raw unnormalized data (default =TRUE)}

\item{norm}{normalize the seurat data (default = TRUE)}
}
\description{
convert BioData to Seurat object but keeping the BioData normalized information.
}
