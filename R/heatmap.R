#' @name heatmap
#' @aliases heatmap,tRNAMINT-method
#' @rdname heatmap-methods
#' @docType methods
#' @description heatmap specific for tRNAMINT data. These heatmaps should be run on size ordered objects 
#' in order to visualize the read counts of the different tRNA fragments starting from the 5' over the center to the 3' and full length fragments.
#' @param x the tRNAMINT object
#' @param colGroup what do you want to have coloured at the top 
#' @param norm.type any of 'all reads' 'all tRNA reads'  or 'Unnormalized'
#' @param codon which codon to focus on (a column in the annotation table)
#' @param fname the output filename
#' @param main the title string for the heatmap default="Heatmap"
#' @title description of function heatmap
#' @export 
setGeneric('heatmap', ## Name
	function ( x, colGroup, norm.type=NULL, codon=NULL, fname=NULL, main="Heatmap") { ## Argumente der generischen Funktion
		standardGeneric('heatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('heatmap', signature = c ('tRNAMINT'),
	definition = function ( x, colGroup, norm.type=NULL, codon=NULL, fname=NULL, main="Heatmap") {
	colors_4(x,colGroup) # just make sure they are defined globaly
	colors_4(x,'tRF.type.s.') # just make sure they are defined globaly
	x <- x$clone() ## I might drop a lot...
	if ( !is.null(norm.type) ) {
		reduceTo(x,what='col', colnames(x$data)[ grep( norm.type, x$samples$NormalizationMode ) ] )
	}
	if ( ! is.null(codon) ) {
		reduceTo(x,what='row', rownames(x$data)[ where( x$annotation[,codon] == 1 ) ] )
	}
	complexHeatmap(x, colGroups=colGroup, ofile=fname,rowGroups=c('tRF.type.s.'), pdf=T )
	
} )
