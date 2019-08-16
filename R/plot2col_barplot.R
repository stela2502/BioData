#' @name plot2col_barplot
#' @aliases plot2col_barplot,BioData-method
#' @rdname plot2col_barplot-methods
#' @docType methods
#' @description Plot n1 vs n2 columns as colored barplot to show the distribution of n2 inside of n1
#' @param x the BioData object
#' @param n1 the samples column for the x axis
#' @param n2 the samples column for the y axis and the coloring of the bars#' @name plot2col_barplot
#' @param file the outfile for this plot default=NULL
#' @param svg create a svg file default=F (png)
#' @param pdf create a pdf default=F (png)
#' @param xScaleFactor rescale width for very long factor names take higher value ( default=10 )
#' @param X11type a type for the plot devices default='cairo'
#' @title description of function plot2col_barplot
#' @export 
if ( ! isGeneric('plot2col_barplot') ){setGeneric('plot2col_barplot', ## Name
	function ( x, n1, n2, file=NULL, svg=F, pdf=F, X11type='cairo', xScaleFactor=10 ) { 
		standardGeneric('plot2col_barplot')
	}
) }

setMethod('plot2col_barplot', signature = c ('BioData'),
	definition = function ( x, n1, n2, file=NULL, svg=F, pdf=F, X11type='cairo', xScaleFactor=10  ) {
	
	data = t((
		table( 
			x$samples[, n1], x$samples[, n2]) /
			as.vector( apply( table( x$samples[, n1], x$samples[, n2]),1, sum))
			*100
			))
	colnames(data) = paste(colnames(data),"\n(n=", sep="")
	colnames(data) = paste(colnames(data), table( x$samples[, n1] ), ")", sep="")
	if ( ! is.null(file) ) {
		file = file.path(x$outpath, paste( collapse='_',
						unlist(strsplit( c(file, n1,'_vs_', n2), '\\s+', perl=T))))
		w = xScaleFactor * ceiling( ncol(data) / 6 )
		if ( svg ) {
			RSvgDevice::devSVG( file=paste(file,'svg',sep='.'), width= w, height=xScaleFactor )
		}
		else if ( pdf ) {
			grDevices::pdf( file=paste(file, 'pdf', sep='.'), width= w, height= xScaleFactor)
		}
		else {
			grDevices::png(file=paste(file, 'png', sep='.'), width= w*100, height=xScaleFactor * 100, type=X11type )
		}
	}
	graphics::barplot( data , col= x$usedObj$colorRange[[n2]], 
			ylab=paste( n2,"[%]"), main=paste( n1, 'vs', n2) )
	graphics::legend("topleft", rev(levels(x$samples[, n2])), 
			fill=rev(x$usedObj$colorRange[[n2]]),bg='white')
	if ( ! is.null(file) ) {
		grDevices::dev.off()
	}
	invisible(x)
} )
