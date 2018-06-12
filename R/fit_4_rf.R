#' @name fit_4_rf
#' @aliases fit_4_rf,BioData-method
#' @rdname fit_4_rf-methods
#' @docType methods
#' @description fit_4_rf removes raw and zscored data and all -1 values from the dat object
#' @param x The BioData object
#' @param copy Return a copy of the object default= TRUE
#' @title Reduce memory footprint and -1 values from a BioData object
#' @export 
setGeneric('fit_4_rf', ## Name
	function ( x, copy= TRUE ) { ## Argumente der generischen Funktion
		standardGeneric('fit_4_rf') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('fit_4_rf', signature = c ('BioData'),
	definition = function ( x, copy= TRUE ) {
	if ( copy ) {
		x <- x$clone()
	}
	x$raw = NULL
	x$zscored = NULL
	m <- min(x$dat)
	if ( m == -1 ) {
		lapply( 1:ncol(x$dat), function(i) { 
					to0 <- which(x$dat[,i] == -1 ) 
					if ( length(to0) > 0 ){
						x$dat[to0,i] = 0
					}
					NULL
					 } 
			 )
	}
	invisible(x)
} )
