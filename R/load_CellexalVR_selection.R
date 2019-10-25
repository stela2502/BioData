#' @name load_CellexalVR_selection
#' @aliases load_CellexalVR_selection,BioData-method
#' @rdname load_CellexalVR_selection-methods
#' @docType methods
#' @description load a selection.txt file created during a CellexalVR run
#' @param x the BioData object
#' @param file the cellexal selection<id>.txt file
#' @param gname the sample group name this selection should get in the BioData object
#' @title description of function load_CellexalVR_selection
#' @export 
if ( ! isGeneric('load_CellexalVR_selection') ){setGeneric('load_CellexalVR_selection', ## Name
	function ( x, file, gname ) { 
		standardGeneric('load_CellexalVR_selection')
	}
) }

setMethod('load_CellexalVR_selection', signature = c ('BioData'),
	definition = function ( x, file, gname ) {
	if ( ! file.exists( file )) {
		stop(paste("Sorry I can not open the file", file) )
	}
	info = read.delim( file, header=F )
	x$samples[,gname] = NA
	OK = match( info[,1], colnames(x) )
	x$samples[OK,gname] = info[,4] +1
	x$usedObj$colorRange[[gname]] = unique(  info[,2] )
	invisible(x)
} )
