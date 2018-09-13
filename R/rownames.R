#' @name rownames
#' @aliases rownames,BioData-method
#' @rdname rownames-methods
#' @docType methods
#' @description wrapper for rownames(x$dat)
#' @param x  TEXT MISSING
#' @param do.NULL see base::colnames (default= TRUE)
#' @param prefix  see base::colnames  T default= "row"
#' @title wrapper for rownames(x$dat)
#' @export 
if ( ! isGeneric('rownames') ){ setGeneric('rownames', ## Name
	function ( x, do.NULL = TRUE, prefix = "row") { ## Argumente der generischen Funktion
		standardGeneric('rownames') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)}

setMethod('rownames', signature = c ('BioData'),
	definition = function ( x, do.NULL = TRUE, prefix = "row") {
	rownames(x$dat, do.NULL = do.NULL, prefix = prefix)
} )
