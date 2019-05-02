getTopSigGenes <- function(x, name, n=100 ) {
	
	all_markers = x$stats[[name]]
	
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
}