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
	function ( x, what='row', colname ) { 
		standardGeneric('changeNames')
	}
)

setMethod('changeNames', signature = c ('BioData'),
	definition = function ( x, what='row', colname ) {
		rep_missing <- function ( v ) {
			notOK <- which(is.na(v))
			if ( length(notOK) > 0) {
				v[notOK] <- sprintf( paste("No--%0",nchar(max(notOK)),"d",sep=""),1:length(notOK))
			}
			v
		}
	if ( what=='row' ) {
		id = match(colname, colnames(x$annotation) )
		if ( is.na(id) ){
			stop( paste( "The colname", colname,"can not be found in the annotation table") )
		}
		if ( length(which(table(as.character(x$annotation[,colname] ) ) > 1 ) ) > 0 ){
			x$rownamescol <- paste(colname, 'unique', sep='.' )
			rownames(x$annotation) <- x$annotation[ , x$rownamescol ] <- x$forceAbsoluteUniqueSample( rep_missing(as.character(x$annotation[,colname] )) )
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
			rownames(x$samples) <- x$samples[ , x$sampleNamesCol ] <- x$forceAbsoluteUniquerep_missing(Sample( as.character(x$samples[,colname] )) )
			
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
