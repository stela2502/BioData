getNgenesFromMultiGroupStats <- function( all_markers, num.sig=1000 ) {

	genes_list <- split( as.vector(all_markers[,'gene']), all_markers[,'cluster'] )
	ret_genes =  ceiling(num.sig / length(table(grp.vec)))
		
	if ( ret_genes < 1)
		ret_genes = 1
	
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
		
	## likely not the best approach..
	deg.genes = NULL
	ret_genes = ret_genes -1
	i = 0
	while ( length( deg.genes ) < num.sig ) {
		ret_genes = ret_genes +1
		i = i+1
		deg.genes = unique(unlist( lapply( genes_list,top_genes ) ))
		bad = which(is.na(deg.genes))
		if ( length(bad) > 0) 
			deg.genes = deg.genes[-bad]
		if ( i > 20)
			break
	}
	
	deg.genes = rownames(x@data)[ match( make.names(deg.genes), make.names( rownames( x@data) ) )]
	deg.genes
}