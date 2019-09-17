#' @name removeNeg
#' @aliases removeNeg,BioData-method
#' @rdname removeNeg-methods
#' @docType methods
#' @description removes -1 values from a SingleCells object
#' @param x teh SingleCells object
#' @param from which slot to remove the data from (dat or zscored) default='dat'
#' @title description of function removeNeg
#' @export 
if ( ! isGeneric('removeNeg') ){setGeneric('removeNeg', ## Name
	function ( x, from='dat' ) { 
		standardGeneric('removeNeg')
	}
) }

setMethod('removeNeg', signature = c ('BioData'),
	definition = function ( x, from='dat' ) {
	if ( from == 'dat'){
		bad = which(x$dat@x == -1)
		if ( length(bad) >0 ){
			x$dat@x[bad]=0
			drop0(x$dat)
		}
	}else if (  from == 'zscored'){
		bad = which(x$zscored@x == -1)
		if ( length(bad) >0 ){
			x$zscored@x[bad]=0
			drop0(x$zscored)
		}
	}
	else {
		stop( paste( 'from', from, 'is not implemented :-('))
	}
	invisible(x)
} )
