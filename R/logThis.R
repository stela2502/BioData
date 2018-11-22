#' @name logThis
#' @aliases logThis,BioData-method
#' @rdname logThis-methods
#' @docType methods
#' @description calculate the log for all data +1
#' @param x the biodata object
#' @title description of function log
#' @export 
if ( ! isGeneric('logThis') ){ setGeneric('logThis', ## Name
	function (x) { 
		standardGeneric('logThis')
	}
)
}else {
	print ("Onload warn generic function 'logThis' already defined - no overloading here!")
}

setMethod('logThis', signature = c ('BioData'),
	definition = function (x) {
	if ( ! x$logged ) {
		#pb <- progress_estimated(100)
		#steps = ceiling(ncol(x$dat)/100)
		if ( is.null(x$raw) ) {
			x$raw = x$dat
		}
		OK = which(x$dat@x > 0)
		x$dat@x[OK] = log(x$dat@x[OK])
		x$logged = TRUE
		
	}
	invisible(x)
} )
