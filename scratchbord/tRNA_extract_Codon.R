heatmap <- function( x, colGroup, norm.type=NULL, codon=NULL, fname=NULL, main="Heatmap") {
	colors_4(x,colGroup) # just make sure they are defined globaly
	x <- x$clone() ## I might drop a lot...
	if ( !is.null(norm.type) ) {
		reduceTo(x,where='col', colnames(x$data)[ grep( norm.type, x$samples$NormalizationMode ) ] )
	}
	if ( ! is.null(codon) ) {
		reduceTo(x,where='row', rownames(x$data)[ where( x$annotation[,codon] == 1 ) ] )
	}
	complexHeatmap(x, colGroups=colGroup, ofile=fname,rowGroups=c('tRF.type.s.'), pdf=T )
	
}