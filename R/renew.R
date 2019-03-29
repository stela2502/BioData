#' @name renew
#' @aliases renew,BioData-method
#' @rdname renew-methods
#' @docType methods
#' @description update the class definition by re-creating the instance
#' @param x the object you want to update
#' @title description of function renew
#' @export 
if ( ! isGeneric('renew') ){ methods::setGeneric('renew', ## Name
	function ( x ) { 
		standardGeneric('renew')
	}
)
}else {
	print ("Onload warn generic function 'renew' already defined - no overloading here!")
}

setMethod('renew', signature = c ('BioData'),
	definition = function ( x ) {
		oldClass <- class(x)
		if ( is.function(x$data) ){
			ret <- BioData$new( 
					dat = x$dat, 
					annotation = x$annotation, 
					Samples = x$samples, 
					name=x$name, 
					namecol=x$sampleNamesCol, 
					namerow= x$rownamescol, 
					outpath=x$outpath 
			)
		}else if ( x$version != utils::sessionInfo('BioData')$otherPkgs$BioData$Version ){
			ret <- BioData$new( 
				dat = x$dat, 
				annotation= x$annotation,
				Samples = x$samples, 
				name=x$name, 
				namecol=x$sampleNamesCol, 
				namerow= x$rownamescol, 
				outpath=x$outpath 
			)
		}
		ret$raw <- x$raw
	ret$usedObj <- x$usedObj
	ret$stats <- x$stats
	ret$snorm <- x$snorm
	ret$zscored <- x$zscored
	x <- ret
	class(x) <- oldClass
	invisible(x)
} )
