% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapse.R
\docType{methods}
\name{collapse}
\alias{collapse}
\alias{collapse,BioData-method}
\title{description of function collapse}
\usage{
collapse(
  x,
  what = "row",
  group,
  fun = function(x) {     mean(x, na.rm = TRUE) }
)
}
\arguments{
\item{x}{the BioData object}

\item{what}{collapse on a row or column data default='row'}

\item{group}{the colnames of the annotaion or samples table}

\item{fun}{the collapsing function (default function(x) {mean(x, ns.rm=TRUE )} )}
}
\description{
This function collapses a dataset on any row or column information using the function fun.
}
