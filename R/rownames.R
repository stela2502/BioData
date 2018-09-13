#' @name rownames
#' @aliases rownames,BioData-method
#' @rdname rownames-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param do.NULL  TEXT MISSING default= TRUE
#' @param prefix  TEXT MISSING default= "row"
#' @title description of function rownames
#' @export 
setGeneric('rownames', ## Name
	function ( x, do.NULL = TRUE, prefix = "row") { ## Argumente der generischen Funktion
		standardGeneric('rownames') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('rownames', signature = c ('BioData'),
	definition = function ( x, do.NULL = TRUE, prefix = "row") {
	rownames(x$dat)
} )
