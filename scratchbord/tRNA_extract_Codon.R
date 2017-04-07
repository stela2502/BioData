plot.codonUsage <- function( x, sname, codons=NULL, min_reads=1, tRF.reliability=NULL, tRF.type=NULL ,main="Pie Chart", color=NULL, fname=NULL ){
	data <- sampleCodonUsage(x, sname, codons=codons, min_reads=min_reads, tRF.reliability=tRF.reliability,tRF.type=tRF.type  )
	if (!is.null(fname)) {
		png(file=paste(x$outpath,fname,"_pie.png", sep=''), width=1000, height=1000)
		pie( data, labels=codons, main=main, col=color)
		dev.off()
		png(file=paste(x$outpath,fname,"_bars.png", sep=''), width=1000, height=1000)
		par(las=2)
		barplot( data, names.arg=codons, main=main, col=color)
		dev.off()
	}
	else {
		pie( data, labels=codons, main=main, col=color)
		x11()
		par(las=2)
		barplot( data, names.arg=codons, main=main, col=color)
	}
	invisible(x)
}