#' @name tRF.typeHeatmap
#' @aliases tRF.typeHeatmap,tRNAMINT-method
#' @rdname tRF.typeHeatmap-methods
#' @docType methods
#' @description tRF.typeHeatmap specific for tRNAMINT data. These heatmaps should be run on size ordered objects 
#' in order to visualize the read counts of the different tRNA fragments starting from the 5' over the center to the 3' and full length fragments.
#' @param x the tRNAMINT object
#' @param colGroup what do you want to have coloured at the top 
#' @param norm.type any of 'all reads' 'all tRNA reads'  or 'Unnormalized'
#' @param tRF.type which tRF.type to focus on (a column in the annotation table)
#' @param fname the output filename
#' @param main the title string for the heatmap default="Heatmap"
#' @param fun collapse the data before plotting; e.g.
#' @param z.score z.score the data after the summary (default = F)
#' @param brks how many shades of heatmap colors ( default 10 )
#' function(x) \{ collaps(x,what='row',group='frag.type.and.length', fun = function(x) \{ sum( x, na.rm=TRUE) \} ) \}  
#' @title description of function heatmap
#' @export 
if ( ! isGeneric('tRF.typeHeatmap') ){ methods::setGeneric('tRF.typeHeatmap', ## Name
	function ( x, colGroup, norm.type=NULL, tRF.type=NULL, fname=NULL, main="Heatmap", fun=function(x) {collapse2codons(x)}, z.score=FALSE, brks=10) { 
		standardGeneric('tRF.typeHeatmap')
	}
)
}else {
	print ("Onload warn generic function 'tRF.typeHeatmap' already defined - no overloading here!")
}

setMethod('tRF.typeHeatmap', signature = c ('tRNAMINT'),
	definition = function ( x, colGroup, norm.type=NULL, tRF.type=NULL, fname=NULL, main="Heatmap", fun=function(x) {collapse2codons(x)}, z.score=FALSE, brks=10) {
		rowGroup=c('Codon')
		colors_4(x,colGroup) # just make sure they are defined globaly
	colors_4(x,'tRF.type.s.') # just make sure they are defined globaly
	x <- x$clone() ## I might drop a lot...
	if ( !is.null(norm.type) ) {
		reduceTo(x,what='col', colnames(x$dat)[ grep( norm.type, x$samples$NormalizationMode ) ] )
	}
	if ( ! is.null(tRF.type) ) {
		reduceTo(x,what='row', rownames(x$dat)[ which( x$annotation$tRF.type.s. == tRF.type ) ] )
	}
	if ( ! is.null(fun)){
		x = fun(x)
	}
	if ( z.score ) {
		z.score(x)
	}
	colors_4(x,colGroup)
	complexHeatmap(x, colGroups=colGroup, ofile=fname,rowGroups=rowGroup, pdf=T, main=main, brks=brks )
	x
} )
