#' @name changeNames
#' @aliases changeNames,BioData-method
#' @rdname changeNames-methods
#' @docType methods
#' @description rename row- or colnames of the data tables
#' @param x the BioData object
#' @param what which dimension names to change ('row' or 'col'; default='row')
#' @param colname a colname in the annotaion or samples table
#' @title description of function changeNames
#' @export 
setGeneric('changeNames', ## Name
	function ( x, what='row', colname ) { ## Argumente der generischen Funktion
		standardGeneric('changeNames') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('changeNames', signature = c ('BioData'),
	definition = function ( x, what='row', colname ) {
	if ( what=='row' ) {
		id = match(colname, colnames(x$annotation) )
		if ( is.na(id) ){
			stop( paste( "The colname", colname,"can not be found in the annotation table") )
		}
		if ( length(which(table(as.character(x$annotation[,colname] ) ) > 1 ) ) > 0 ){
			x$rownamescol <- paste(colname, 'unique', sep='.' )
			x$annotation[ , x$rownamescol ] <- x$forceAbsoluteUniqueSample( as.character(x$annotation[,colname] ) )
		}else {
			x$rownamescol = colname
		}
		rownames(x$dat) <- as.character(x$annotation[ , x$rownamescol ]) 
		if ( ! is.null(x$raw) ) {
			rownames(x$raw) <- as.character(x$annotation[ , x$rownamescol ]) 
		}
		if ( ! is.null(x$zscored) ) {
			rownames(x$zscored) <- as.character(x$annotation[ , x$rownamescol ]) 
		}
	}
	else if ( what == 'col') {
		id = match(colname, colnames(x$samples) )
		if ( is.na(id) ){
			stop( paste( "The colname", colname,"can not be found in the samples table") )
		}
		if ( length(which(table(as.character(x$samples[,colname] ) ) > 1 ) ) > 0 ){
			x$sampleNamesCol <- paste(colname, 'unique', sep='.' )
			x$samples[ , x$sampleNamesCol ] <- x$forceAbsoluteUniqueSample( as.character(x$samples[,colname] ) )
		}else {
			x$sampleNamesCol = colname
		}
		colnames(x$dat) <- as.character(x$samples[ , x$sampleNamesCol ]) 
		if ( ! is.null(x$raw) ) {
			colnames(x$raw) <- as.character(x$samples[ , x$sampleNamesCol ]) 
		}
		if ( ! is.null(x$zscored) ) {
			colnames(x$zscored) <- as.character(x$samples[ , x$sampleNamesCol ]) 
		}
	}
	invisible(x)
} )
