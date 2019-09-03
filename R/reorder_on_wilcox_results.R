#' @name reorder_on_wilcox_results
#' @aliases reorder_on_wilcox_results,BioData-method
#' @rdname reorder_on_wilcox_results-methods
#' @docType methods
#' @description uses a wilcox result table gene overlap between the groups to identify close groups.
#' @param x the BioData object
#' @param statsN stats table to use
#' @param column samples column to change
#' @param cutoff p value cut of for the wilcox tables default=1e-8
#' @title description of function reorder_on_wilcox_results
#' @export 
if ( ! isGeneric('reorder_on_wilcox_results') ){setGeneric('reorder_on_wilcox_results', ## Name
	function ( x, statsN, column , cutoff=1e-8 ) { 
		standardGeneric('reorder_on_wilcox_results')
	}
) }

setMethod('reorder_on_wilcox_results', signature = c ('BioData'),
	definition = function ( x, statsN, column , cutoff=1e-8 ) {
	if ( is.null(x$stats[[statsN]])) {
		stop( paste( seq="", "statsN '",statsN,"' is not in the list of stats results:", 
						paste(collapse="', '", names(x$stats))))
	}
	if ( is.null( x$samples[,column])) {
		stop( paste ( "Samples column", column, "is not defined"))
	}
	distanceM = NULL
	data = x$stats[[statsN]]
	data = data[which(data[,'p_val_adj'] < cutoff),]
	for ( k in names(table(data[,'cluster'])) ){
		this = as.vector(data[which(data[,'cluster'] == k),'gene'])
		
		distanceM = rbind(distanceM, unlist( lapply( names(table(data[,'cluster'])) , function(i, other) {
									this = as.vector(data[which(data[,'cluster'] == i),'gene'])
									length( intersect( this, other))
								}, this ) ) / length(this)
		)
	}
	diag(distanceM) = 0
	rownames(distanceM) = colnames(distanceM) = levels(x$samples[,column])
	rownames(distanceM)
	new_order =  levels(x$samples[,column])[hclust(dist( distanceM), method='ward.D2' )$order]
	reorder_grouping( x, group = column, new_order= new_order )
	
} )
