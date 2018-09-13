#' @name colnames
#' @aliases colnames,BioData-method
#' @rdname colnames-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param do.NULL  TEXT MISSING default= TRUE
#' @param prefix  TEXT MISSING default= "col"
#' @title description of function colnames
#' @export 
setGeneric('colnames', ## Name
	function ( x, do.NULL = TRUE, prefix = "col") { ## Argumente der generischen Funktion
		standardGeneric('colnames') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('colnames', signature = c ('BioData'),
	definition = function ( x, do.NULL = TRUE, prefix = "col") {
	colnames(x$dat)
} )
