#' @name add_to_stat
#' @aliases add_to_stat,BioData-method
#' @rdname add_to_stat-methods
#' @docType methods
#' @description Adds a stst table to the object
#' @param x The BioDat object
#' @param stat the stats table
#' @param name the name for this stats table
#' @title Add a stsitics table to the $stats slot
#' @export 
setGeneric('add_to_stat', ## Name
	function ( x, stat, name ) { 
		standardGeneric('add_to_stat')
	}
)

setMethod('add_to_stat', signature = c ('BioData'),
	definition = function ( x, stat, name ) {
	if ( length( match('p.adj BF',colnames(stat) )) == 0 & length( match('hurdle',colnames(stat) ))) {
		stat = cbind( stat, 'p.adj BF' = p.adjust(stat[,'hurdle'], method='bonferroni') )
	}
	if ( ! is.na( match( name, names(x$stats)))){
		x$stats[[ match( name, names(x$stats)) ]] <- stat
	}else {
		x$stats[[ length( x$stats ) +1 ]] <- stat
		names(x$stats)[length(x$stats) ] <- name
	}
	invisible(x)
} )
