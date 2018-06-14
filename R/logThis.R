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
		if ( is.null(x$raw) ) {
			x$raw = x$dat
		}
		rn = rownames(x$dat)
		cn = colnames(x$dat)
		x$dat <- Matrix(apply( x$dat, 2, function(x){
					ok = which( x > -1 )
					x[ok] = log(x[ok]+1)
					x
				} ))
		rownames(x$dat) = rn
		colnames(x$dat) = cn
		x$logged = TRUE
	}
	invisible(x)
} )
