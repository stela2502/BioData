#' @name clusters_gene
#' @aliases clusters_gene,BioData-method
#' @rdname clusters_gene-methods
#' @docType methods
#' @description Clusters the data either based on the raw data values or any MDS data type and adds the grouping into the annotation table.
#' @param dataObj the BioData object
#' @param clusterby cluster on raw data or MDS clustered data default="raw"
#' @param useGrouping do nothing and simply use tis grouping default=NULL
#' @param groups.n how many groups should be detected default= 3
#' @param ctype cluster type - either 'hierarchical clust' or 'kmeans' default = 'hierarchical clust'
#' @param onwhat this option has been kept for the Fluidigm data as there FACS data can also be used default = 'Expression'
#' @param cmethod the method to used with the hclust clustering (default = 'ward.D2')
#' @param name the name for the new grouping (default = 'auto_clusters_gene.1:n')
#' @title description of function clusters_gene
#' @export 
if ( ! isGeneric('clusters_gene') ){ setGeneric('clusters_gene', ## Name
	function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3, 
			ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D2', name=NULL ) {
		standardGeneric('clusters_gene')
	}
)
}else {
	print ("Onload warn generic function 'clusters_gene' already defined - no overloading here!")
}

setMethod('clusters_gene', signature = c ('BioData'),
	definition = function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3,
				ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D2', name=NULL ) {
	
			transpose(dataObj)
			clusters(dataObj, clusterby= clusterby, useGrouping=useGrouping, groups.n = groups.n,
					ctype=ctype, onwhat=onwhat, cmethod=cmethod, name=name)
			transpose(dataObj)
			invisible(dataObj)
} )
