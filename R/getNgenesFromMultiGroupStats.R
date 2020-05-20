#' @name getNgenesFromMultiGroupStats
#' @aliases getNgenesFromMultiGroupStats,BioData-method
#' @rdname getNgenesFromMultiGroupStats-methods
#' @docType methods
#' @description return the top significant genes in an even spacing for all analyzed goups.
#' @param all_markers the data frame containg all results (columns 'gene', 'cluster' are important)
#' @param num.sig how many genes to return (default=1000)
#' @title description of function getNgenesFromMultiGroupStats
#' @export 
if ( ! isGeneric('getNgenesFromMultiGroupStats') ){setGeneric('getNgenesFromMultiGroupStats', ## Name
	function ( all_markers, num.sig=1000 ) { 
		standardGeneric('getNgenesFromMultiGroupStats')
	}
) }

setMethod('getNgenesFromMultiGroupStats', signature = c ('data.frame'),
	definition = function ( all_markers, num.sig=1000 ) {

	grp.vec =unique( all_markers[,'cluster'])
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
	m = match( all_markers[,'gene'], deg.genes)
	t = table( all_markers[which(!is.na(m)),'cluster'])
	v = sd(t)
	m = mean(t)

	print(paste("selected mean ~",round(m,2), "genes from the", length(grp.vec), 
		"groups with a sd of ~", round(v,2), "(", round(1- length(deg.genes) /(m*length(grp.vec)),2),"% overlapping)"))
	#deg.genes = 
	deg.genes
} )
