#' @name colnames
#' @aliases colnames,BioData-method
#' @rdname colnames-methods
#' @docType methods
#' @description simple wrapper around colnames(x$dat)
#' @param x the BioData object
#' @param do.NULL see base::colnames (default= TRUE)
#' @param prefix see base::colnames  T default= "col"
#' @title warpper for colnames(x$dat)
#' @export 
if ( ! isGeneric('colnames') ){setGeneric('colnames', ## Name
	function ( x, do.NULL = TRUE, prefix = "col") { ## Argumente der generischen Funktion
		standardGeneric('colnames') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)}

setMethod('colnames', signature = c ('BioData'),
	definition = function ( x, do.NULL = TRUE, prefix = "col") {
	colnames(x$dat, do.NULL = do.NULL, prefix = prefix)
} )
