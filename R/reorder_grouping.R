#' @name reorder_grouping
#' @aliases reorder_grouping,BioData-method
#' @rdname reorder_grouping-methods
#' @docType methods
#' @description Reorder a grouping after some checks that the id string is suitable to do this
#' @param obj the BioDat object
#' @param group the sample group to re-order
#' @param new_order the new order (ordered factor levels)
#' @param what samples or annotation column to re-order? (default 'col')
#' @title description of function reorder_grouping
#' @export 
setGeneric('reorder_grouping', ## Name
	function (obj, group, new_order, what='col' ) { ## Argumente der generischen Funktion
		standardGeneric('reorder_grouping') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('reorder_grouping', signature = c ('BioData'),
	definition = function (obj, group, new_order, what='col' ) {
		if ( what == 'col '){
			k <- length(levels(obj$samples[,group]) )
			err = NULL
			if ( ! all.equal( sort(as.character(new_order)), sort(levels(obj$samples[,group])) ) == T ) {
				stop( paste("Missing/duplicated levels!", paste(sort(new_order), collapse=", "), sep="\n"))
			}
			obj$samples[,group] = factor( as.vector(obj$samples[,group]), levels=new_order )
			print ("sample group order changed")
		}else {
			k <- length(levels(obj$annotation[,group]) )
			err = NULL
			if ( ! all.equal( sort(as.character(new_order)), sort(levels(obj$annotation[,group])) ) == T ) {
				stop( paste("Missing/duplicated levels!", paste(sort(new_order), collapse=", "), sep="\n"))
			}
			obj$annotation[,group] = factor( as.vector(obj$annotation[,group]), levels=new_order )
			print ("annotation group order changed")
		}
		invisible(obj)
	} 
)
