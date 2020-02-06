#' @name plotCodonUsage
#' @aliases plotCodonUsage,BioData-method
#' @rdname plotCodonUsage-methods
#' @docType methods
#' @description plots the codon usage for one sample
#' The function uses sampleCodonUsage() to extract this information from for the sample
#' @param x the tRNAMINT object
#' @param sname the sample you want to depict
#' @param codons the codon list to use default=all
#' @param min_reads use all tRNA fragments with at least that many reads for the sample default=1
#' @param tRF.reliability MINT reports an exclusive and ambiguous set of tRNA fragments state one of them or use both default=NULL == both
#' @param tRF.type MINT reports different types of tRNA fragments (3'-half, 5'-half, ...) default=NULL == use all
#' @param main the figure title default="Pie Chart"
#' @param color a vector of color information to use for the codons default= rainbow
#' @param fname the filename to plot to as png (default use x11 as plotting device)
#' @title description of function plotCodonUsage
#' @export 
if ( ! isGeneric('plotCodonUsage') ){ methods::setGeneric('plotCodonUsage', ## Name
	function ( x, sname, codons=NULL, min_reads=1, tRF.reliability=NULL, tRF.type=NULL ,main="Pie Chart", color=NULL, fname=NULL ) { 
		standardGeneric('plotCodonUsage')
	}
)
}else {
	print ("Onload warn generic function 'plotCodonUsage' already defined - no overloading here!")
}

setMethod('plotCodonUsage', signature = c ('BioData'),
	definition = function ( x, sname, codons=NULL, min_reads=1, tRF.reliability=NULL, tRF.type=NULL ,main="Pie Chart", color=NULL, fname=NULL ) {
	data <- sampleCodonUsage(x, sname, codons=codons, min_reads=min_reads, tRF.reliability=tRF.reliability,tRF.type=tRF.type  )
	if ( is.null(color) ) {
		color = grDevices::rainbow(length(data))
	}
	if (!is.null(fname)) {
		grDevices::png(file=paste(x$outpath,fname,"_pie.png", sep=''), width=1000, height=1000)
		graphics::pie( data, labels=codons, main=main, col=color)
		grDevices::dev.off()
		grDevices::png(file=paste(x$outpath,fname,"_bars.png", sep=''), width=1000, height=1000)
		graphics::par(las=2)
		graphics::barplot( data, names.arg=codons, main=main, col=color)
		grDevices::dev.off()
	}
	else {
		graphics::pie( data, labels=codons, main=main, col=color)
		grDevices::x11()
		graphics::par(las=2)
		graphics::barplot( data, names.arg=codons, main=main, col=color)
	}
	invisible(x)
} )
