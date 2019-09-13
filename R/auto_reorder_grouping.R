#' @name auto_reorder_grouping
#' @aliases auto_reorder_grouping,BioData-method
#' @rdname auto_reorder_grouping-methods
#' @docType methods
#' @description reorder the grouping of single cell data creating assumed summary data
#' @param x the SingleCells object
#' @param group The grouping you wat to reorder
#' @title automaticly reorder a grouping based on summmary expression hclust
#' @export 
if ( ! isGeneric('auto_reorder_grouping') ){setGeneric('auto_reorder_grouping', ## Name
	function ( x, group ) { 
		standardGeneric('auto_reorder_grouping')
	}
) }

setMethod('auto_reorder_grouping', signature = c ('SingleCells'),
	definition = function ( x, group ) {
	
	colapsed = collaps( x, groupCol= group , copy=T, by='sum')
	class(colapsed) = c('BioData', 'R6')
	colapsed$raw= colapsed$dat
	normalize(colapsed, force=T)
	
	clusters( colapsed, groups.n= 2 )
	new_order <- colapsed$usedObj$hc$labels[colapsed$usedObj$hc$order]
	reorder_grouping( x, group= group, new_order= new_order )
	new_order
} )
