#' @name getTopSigGenes
#' @aliases getTopSigGenes,BioData-method
#' @rdname getTopSigGenes-methods
#' @docType methods
#' @description obtain a list of most significant genes per comparison froma Seurat or Cpp stats test.
#' @param x the BioData object
#' @param name the grouping name (samples column name)
#' @param n the amount of genes per stat test patrt (default 100)
#' @title obtain a list of most significant genes per comparison froma Seurat or Cpp stats test.
#' @export 
setGeneric('getTopSigGenes', ## Name
	function (x, name, n=100 ) { ## Argumente der generischen Funktion
		standardGeneric('getTopSigGenes') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('getTopSigGenes', signature = c ('BioData'),
	definition = function (x, name, n=100 ) {
	
	all_markers = x$stats[[name]]
	all_markers = all_markers[ which( all_markers$logFC > 0),]
	
	if ( length( match( 'gene', colnames(all_markers))) == 0  ){
		stop("This function is mean to parse Cpp or Seurat stat output tables only" )
	}
	
	genes_list <- split( as.vector(all_markers[,'gene']), all_markers[,'cluster'] )
	ret_genes = n
	top_genes <- function( x ) {
		if ( length(x) == 0) {
			NA
		}
		else if ( length(x) < ret_genes ) {
			x
		}else {
			x[1:ret_genes]
		}
	}
	
	deg.genes = NULL
	deg.genes = unique(unlist( lapply( genes_list,top_genes ) ))
	bad = which(is.na(deg.genes))
	if ( length(bad) > 0)
		deg.genes = deg.genes[-bad]
	deg.genes
} )
