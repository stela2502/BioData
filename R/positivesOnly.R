#' @name positivesOnly
#' @aliases positivesOnly,BioData-method
#' @rdname positivesOnly-methods
#' @docType methods
#' @description get a stored Matrix data set without the normalization drop out information
#' @param x the SingleCells object
#' @param what which slot to get default='dat'
#' @title description of function positivesOnly
#' @export 
if ( ! isGeneric('positivesOnly') ){setGeneric('positivesOnly', ## Name
	function (x, what='dat' ) { 
		standardGeneric('positivesOnly')
	}
) }

setMethod('positivesOnly', signature = c ('BioData'),
	definition = function (x, what='dat' ) {
	if(  is.null( x[[what]]) ) {
		stop(paste("data type", what, "is not defined"))
	}
	ret = x[[what]]
	if (what== 'zscored') {
		bad = which(ret@x == -20 )
	}else {
		bad = which(ret@x == -1 )
	}
	if ( length(bad) > 0 ) {
		ret@x[bad] = 0
		ret = drop0(ret)
	}
	ret
} )
