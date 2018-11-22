#' @name CodonHeatmap
#' @aliases CodonHeatmap,tRNAMINT-method
#' @rdname CodonHeatmap-methods
#' @docType methods
#' @description CodonHeatmap specific for tRNAMINT data. These heatmaps should be run on size ordered objects 
#' in order to visualize the read counts of the different tRNA fragments starting from the 5' over the center to the 3' and full length fragments.
#' @param x the tRNAMINT object
#' @param colGroup what do you want to have coloured at the top 
#' @param rowGroup a list of annotaon table columns (default =  'tRF.type.s.')
#' @param norm.type any of 'all reads' 'all tRNA reads'  or 'Unnormalized'
#' @param codon which codon to focus on (a column in the annotation table)
#' @param fname the output filename
#' @param main the title string for the heatmap default="Heatmap"
#' @param fun collapse the data before plotting; e.g.
#' function(x) \{ collaps(x,what='row',group='frag.type.and.length', fun = function(x) \{ sum( x, na.rm=TRUE) \} ) \}  
#' @title description of function heatmap
#' @export 
if ( ! isGeneric('CodonHeatmap') ){ setGeneric('CodonHeatmap', ## Name
	function ( x, colGroup,  rowGroup=c('tRF.type.s.'), norm.type=NULL, codon=NULL, fname=NULL, main="Heatmap", fun=NULL) { 
		standardGeneric('CodonHeatmap')
	}
)
}else {
	print ("Onload warn generic function 'CodonHeatmap' already defined - no overloading here!")
}

setMethod('CodonHeatmap', signature = c ('tRNAMINT'),
	definition = function ( x, colGroup, rowGroup=c('tRF.type.s.'), norm.type=NULL, codon=NULL, fname=NULL, main="Heatmap", fun=NULL) {
	colors_4(x,colGroup) # just make sure they are defined globaly
	colors_4(x,'tRF.type.s.') # just make sure they are defined globaly
	x <- x$clone() ## I might drop a lot...
	if ( !is.null(norm.type) ) {
		reduceTo(x,what='col', colnames(x$dat)[ grep( norm.type, x$samples$NormalizationMode ) ] )
	}
	if ( ! is.null(codon) ) {
		reduceTo(x,what='row', rownames(x$dat)[ which( x$annotation[,codon] == 1 ) ] )
	}
	if ( ! is.null(collapse)){
		fun(x)
	}
	complexHeatmap(x, colGroups=colGroup, ofile=fname,rowGroups=rowGroup, pdf=T )
	
} )
