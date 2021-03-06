#' @name logThis
#' @aliases logThis,BioData-method
#' @rdname logThis-methods
#' @docType methods
#' @description calculate the log for all data +1
#' @param x the biodata object
#' @title description of function log
#' @export 
if ( ! isGeneric('logThis') ){ setGeneric('logThis', ## Name
	function (x) { ## Argumente der generischen Funktion
		standardGeneric('logThis') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'logThis' already defined - no overloading here!")
}

setMethod('logThis', signature = c ('BioData'),
	definition = function (x) {
	if ( ! x$logged ) {
		x$data <- apply( x$data, 2, function(x){ log(x+1) } )
		x$logged = TRUE
	}
	invisible(x)
} )
