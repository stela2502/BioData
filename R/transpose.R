#' @name transpose
#' @aliases transpose,02_transpose-method
#' @rdname transpose-methods
#' @docType methods
#' @description transpose the object DELETING all stats and usedObj!
#' @param x the BioData object
#' @title description of function t
#' @export 
setGeneric('transpose', ## Name
		function (x) { ## Argumente der generischen Funktion
			standardGeneric('transpose') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

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
			x$data <- t(x$data)
			tmp <- x$sampleNamesCol
			x$sampleNamesCol <- x$rownamescol
			x$rownamescol <- tmp
			if ( ! is.null(x$raw) ) {
				x$raw < - t(x$raw)
			}
			tmp = x$usedObj$MDS
			x$usedObj$MDS = x$usedObj$MDSgenes
			x$usedObj$MDSgenes = tmp
			invisible(x)
		} 
)
