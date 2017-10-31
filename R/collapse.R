#' @name collapse
#' @aliases collapse,BioData-method
#' @rdname collapse-methods
#' @docType methods
#' @description This function collapses a dataset on any row or column information using the function fun.
#' @param x the BioData object
#' @param what collapse on a row or column data default='row'
#' @param group the colnames of the annotaion or samples table
#' @param fun the collapsing function (default function(x) {mean(x, ns.rm=TRUE )} )
#' @title description of function collapse
#' @export 
if ( ! isGeneric('collapse') ){ setGeneric('collapse', ## Name
	function ( x, what='row', group, fun = function(x){ mean(x, na.rm=TRUE)} ) { ## Argumente der generischen Funktion
		standardGeneric('collapse') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'collapse' already defined - no overloading here!")
}

setMethod('collapse', signature = c ('BioData'),
	definition = function ( x, what='row', group, fun = function(x){ mean(x, na.rm=TRUE)} ) {
		stop( "Broken function - use collaps instead of collapse!")
		if ( is.null(x$raw) ) {
			x$raw <- x$dat
		}
	if ( what == 'row') {
		if ( is.null( x$annotation[[group]]) ) {
			stop( paste("No colname", group, "defined in the annotation data" ) )
		}
		print ("the annotation data will be corruped as I only use the first row that matches to the grouping data in the new table!" )
		x$usedObj$tmp = data.frame()
		x$dat <- t(sapply( unique( as.character(x$annotation[,group])), function( name ){
			ids <- which(x$annotation[,group] == name)
			x$usedObj$tmp <- rbind(x$usedObj$tmp, x$annotation[ids[1],])
			red <- x$data[ids,]
			apply( red,2,fun )
				} ))
		x$annotation <- x$usedObj$tmp
		x$usedObj$tmp = NULL
		
	}else if ( what=='col') {
		if ( is.null( x$samples[[group]]) ) {
			stop( paste("No colname", group, "defined in the samples data" ) )
		}
		print ("the samples data will be corruped as I only use the first row that matches to the grouping data in the new table!" )
		x$usedObj$tmp = data.frame()
		x$data <- sapply( unique( as.character(x$samples[,group])), function( name ){
							ids <- which(x$samples[,group] == name)
							x$usedObj$tmp <- rbind(x$usedObj$tmp, x$samples[ids[1],])
							red <- x$data[,ids]
							apply( red,1,fun )
						} )
		x$samples <- x$usedObj$tmp
		x$usedObj$tmp = NULL
		
	}else {
		stop( "What has to be either 'row' or 'col'" )
	}
	invisible(x)
} )
