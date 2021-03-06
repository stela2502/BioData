% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusters.R
\docType{methods}
\name{clusters}
\alias{clusters}
\alias{clusters,BioData-method}
\title{create a grouping based on either the raw data or a MDS projection}
\usage{
clusters(
  dataObj,
  clusterby = "raw",
  useGrouping = NULL,
  groups.n = 3,
  ctype = "hierarchical clust",
  onwhat = "Expression",
  cmethod = "ward.D2",
  name = NULL
)
}
\arguments{
\item{dataObj}{the BioData object}

\item{clusterby}{is used to specify the MDS type if onwhat is set to MDS (PCA, TSNE_R ...) or the raw, zscored or raw slot for expression data (default 'raw')}

\item{useGrouping}{do nothing and simply use this grouping default=NULL}

\item{groups.n}{how many groups should be detected default= 3}

\item{ctype}{cluster type - either 'hierarchical clust', 'mclust' or 'kmeans' default = 'hierarchical clust'}

\item{onwhat}{This option selects the source data ('Expression', or 'MDS')}

\item{cmethod}{the method to used with the hclust clustering (default = 'ward.D2')}

\item{name}{the name for the new grouping (default = 'auto_clusters.1:n')}
}
\description{
This function clusters the data based on hclust, mclust or the kmeans algorithm.
This function can be used to cluster any data slot the BioData class offers (see combination of onwhat and clusterby).
}
\examples{
clusters( TestData , clusterby = "TSNE_R", groups.n = 3, ctype = "kmeans", onwhat= 'MDS', name = "kmeansTSNE_R clusters" )
}
