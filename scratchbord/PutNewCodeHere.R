reorder_on_wilcox_results <- function( x, statsN, column , cutoff=1e-8 ) {
	if ( is.null(x$stats[[statsN]])) {
		stop( paste( seq="", "statsN '",statsN,"' is not in the list of stats results:", 
						paste(collapse("', '", names(x$stats)))))
	}
	if ( is.null( x$samples[,column])) {
		stop( paste ( "Samples column", column, "Iis not defined"))
	}
	distanceM = NULL
	data = x$stats[[statsN]]
	for ( k in 1:max(as.numeric(data[,'cluster']))){
		this = as.vector(data[which(data[,'cluster'] == k),'gene'])
		distanceM = rbind(distanceM, unlist( lapply( 1:max(as.numeric(data[,'cluster'])) , function(i, other) {
									this = as.vector(data[which(data[,'cluster'] == i),'gene'])
									length( intersect( this, other))
								}, this ) ) / length(this)
		)
	}
	diag(distanceM) = 0
	new_order =  hclust(dist( distanceM), method='ward.D2' )$order
	reorder_grouping( x, group = column, new_order= levels( x$samples[,column])[new_order])
	
}
