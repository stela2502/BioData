#' @name force.numeric
#' @aliases force.numeric,BioData-method
#' @rdname force.numeric-methods
#' @docType methods
#' @description make sure, that all variables in the data table are numeric. Add NAs if necessary.
#' @param x the BioData object
#' @title description of function force.numeric
#' @export 
if ( ! isGeneric('force.numeric') ){ setGeneric('force.numeric', ## Name
	function ( x ) { ## Argumente der generischen Funktion
		standardGeneric('force.numeric') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'force.numeric' already defined - no overloading here!")
}

setMethod('force.numeric', signature = c ('BioData'),
	definition = function ( x ) {
	for ( i in 1: ncol(x$data) ) {
		if ( ! is.numeric(x$data[,i]) ) {
			x$data[,i] <- as.numeric(as.vector(x$data[,i]))
		}
	}
	invisible(x)
} )
