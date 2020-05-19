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
	function ( x, copy= TRUE ) { 
		standardGeneric('fit_4_rf')
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
		x$dat@x [ which( x$dat@x == -1)] = 0;
		x$dat = Matrix::drop0( x$dat)
	}
	invisible(x)
} )
