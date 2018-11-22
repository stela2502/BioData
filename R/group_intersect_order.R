#' @name group_intersect_order
#' @aliases group_intersect_order,BioData-method
#' @rdname group_intersect_order-methods
#' @docType methods
#' @description use a second ordering to reorder a given sample grouping
#' @param x the BioData object
#' @param group the group to be reordered
#' @param otherGroup the leading grouping
#' @title description of function group_intersect_order
#' @export 
setGeneric('group_intersect_order', ## Name
	function ( x, group, otherGroup ) { 
		standardGeneric('group_intersect_order')
	}
)

setMethod('group_intersect_order', signature = c ('BioData'),
	definition = function ( x, group, otherGroup ) {
	r = lapply( levels(x$samples[, group]), 
			function(n) { 
				t = table(x$samples[, otherGroup][which(x$samples[, group] == n )])
				c(names(t[which( t == max(t))[1]]),n) 
			} 
	)
	r= t(data.frame(r))
	rownames(r) = NULL
	new_order = order(factor( r[,1], levels( x$samples[, otherGroup] )))
	reorder_grouping(x, group= group,  new_order=as.vector(r[new_order,2]),  what='col')
	invisible(x)
}  )
