#' @name transpose
#' @aliases transpose,02_transpose-method
#' @rdname transpose-methods
#' @docType methods
#' @description transpose the object DELETING all stats and usedObj!
#' @param x the BioData object
#' @title description of function t
#' @export 
if ( ! isGeneric('transpose') ){ setGeneric('transpose', ## Name
		function (x) { 
			standardGeneric('transpose')
		}
)
}else {
	print ("Onload warn generic function 'transpose' already defined - no overloading here!")
}

setMethod('transpose', signature = c ('BioData'),
		definition = function (x) {
			if ( is.null(x$usedObj[['transposed']])){
				x$usedObj[['transposed']] = TRUE
			}else {
				x$usedObj[['transposed']] = FALSE
			}
			tmp <- x$annotation
			x$annotation <- x$samples
			x$samples <- tmp
			x$dat <- Matrix::t(x$dat)
			tmp <- x$sampleNamesCol
			x$sampleNamesCol <- x$rownamescol
			x$rownamescol <- tmp
			if ( ! is.null(x$raw) ) {
				x$raw <- Matrix::t(x$raw)
			}
			if ( ! is.null(x$zscored) ) {
				x$zscored <- Matrix::t(x$zscored)
			}
			tmp = x$usedObj$MDS
			x$usedObj$MDS = x$usedObj$MDSgene
			x$usedObj$MDSgene = tmp
			invisible(x)
		} 
)
