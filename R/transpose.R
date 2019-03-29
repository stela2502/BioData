#' @name transpose
#' @aliases transpose,02_transpose-method
#' @rdname transpose-methods
#' @docType methods
#' @description transpose the object DELETING all stats and usedObj!
#' @param x the BioData object
#' @title description of function t
#' @export 
if ( ! isGeneric('transpose') ){ methods::setGeneric('transpose', ## Name
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
			new = list()
			for ( id in grep ( 'MDS', names(x$usedObj)) ) {
				n = names(x$usedObj)[id]
				new_N = NULL
				if ( length(grep('MDSgenes',n)) == 1 ) {
					new_N = stringr::str_replace(n, 'MDSgenes', 'MDS' )
				}else {
					new_N = stringr::str_replace(n, 'MDS', 'MDSgenes' )
				}
				new[[new_N]] = x$usedObj[[n]]
				x$usedObj[[n]] = NULL
			}
			if ( !is.na( match( 'pr', names(x$usedObj)))){
				new[['prGenes']] = x$usedObj[['pr']]
				x$usedObj[['pr']] = NULL
			}
			if ( !is.na( match( 'prGenes', names(x$usedObj)))){
				new[['pr']] = x$usedObj[['prGenes']]
				x$usedObj[['prGenes']] = NULL
			}
			lapply( names(new) , function(n) { x$usedObj[[n]] = new[[n]] })
			rm(new)
		
			invisible(x)
		} 
)
