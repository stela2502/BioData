#' @name plotLegend
#' @aliases plotLegend,BioData-method
#' @rdname plotLegend-methods
#' @docType methods
#' @description Creates the legend for one sample or annotation variable using either the inbuilt 
#' color or an external one in the col option. If a file is given, the data will be plotted into a file. 
#' If neither pdf of svg is true the data will be plotted to a png file. File endings will always be added - so do not do that yourself.
#' @param x the StefansExpressionSet object
#' @param colname the name of a column in the annotation or samples tables that should be described in the legend
#' @param file the optional outfile
#' @param svg create a svg file default=F
#' @param pdf create a pdf file default=F
#' @param col a vector of color names default=NULL
#' @title description of function plot.legend
#' @export 
setGeneric('plotLegend', ## Name
		function ( x, colname, file=NULL, svg=F, pdf=F, col=NULL, X11type='cairo' ) { ## Argumente der generischen Funktion
			standardGeneric('plotLegend') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plotLegend', signature = c ('BioData'),
		definition = function ( x, colname, file=NULL, svg=F, pdf=F, col=NULL, X11type='cairo' ) {
			if ( is.null(col) ){
				col=x$usedObj[['colorRange']][[colname]]
			}
			n =1
			if (! is.na(match(colname, colnames(x$annotation) ))) {
				n <- length( levels( x$annotation[, colname] ) )
			}else if ( ! is.na(match(colname, colnames(x$samples) ))){
				n <- length( levels( x$samples[, colname] ) )
			}else{
				stop ( paste( 'The column name',colname, 
								'is nether defined in the samples nor the annotation table:', 
								paste(colnames(x$samples),collapse=' '),
								paste(x$annotation, collapse=' ')
						))
			}
			if ( ! is.null(file) ) {
				file = file.path(x$outpath, paste( collapse='_',unlist(strsplit( c(file, colname), '\\s+', perl=T))))
				h = 4 * ceiling(n /17)
				if ( svg ) {
					devSVG( file=paste(file,'svg',sep='.'), width= 4, height=h )
				}
				else if ( pdf ) {
					pdf( file=paste(file, 'pdf', sep='.'), width= 4, height=h )
				}
				else {
					png(file=paste(file, 'png', sep='.'), width= 400, height=h*100, type=X11type )
				}
			}
			plot(1, type="n", axes=F, xlab="", ylab="", main=colname)
			if (! is.na(match(colname, colnames(x$annotation) ))) {
				legend( 'top', levels( x$annotation[, colname] ), fill=col )
			}else if ( ! is.na(match(colname, colnames(x$samples) ))){
				legend( 'top', levels( x$samples[, colname] ), fill=col )
			}else{
				stop ( paste( 'The column name',colname, 
								'is nether defined in the samples nor the annotation table:', 
								paste(colnames(x$samples),collapse=' '),
								paste(x$annotation, collapse=' ')
				))
			}
			if ( ! is.null(file) ) {
				dev.off()
			}
		}
)


