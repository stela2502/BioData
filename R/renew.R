#' @name renew
#' @aliases renew,BioData-method
#' @rdname renew-methods
#' @docType methods
#' @description update the class definition by re-creating the instance
#' @param x the object you want to update
#' @title description of function renew
#' @export 
if ( ! isGeneric('renew') ){ setGeneric('renew', ## Name
	function ( x ) { ## Argumente der generischen Funktion
		standardGeneric('renew') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'renew' already defined - no overloading here!")
}

setMethod('renew', signature = c ('BioData'),
	definition = function ( x ) {
	ret <- BioData$new( 
			dat = cbind(x$data, x$annotation), 
			Samples = x$samples, 
			name=x$name, 
			namecol=x$sampleNamesCol, 
			namerow= x$rownamescol, 
			outpath=x$outpath 
	)
	ret$usedObj <- x$usedObj
	ret$stats <- x$stats
	ret$snorm <- x$snorm
	ret$zscored <- x$zscored
	x <- ret
	invisible(x)
} )

setMethod('renew', signature = c ('tRNAMINT'),
		definition = function ( x ) {
			ret <- tRNAMINT$new( 
					dat = cbind(x$data, x$annotation), 
					Samples = x$samples, 
					name=x$name, 
					namecol=x$sampleNamesCol, 
					namerow= x$rownamescol, 
					outpath=x$outpath 
			)
			ret$usedObj <- x$usedObj
			ret$stats <- x$stats
			ret$snorm <- x$snorm
			ret$zscored <- x$zscored
			x <- ret
			invisible(x)
		} )