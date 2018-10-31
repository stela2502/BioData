
CreateBin <- function(x, group, where='sample', colFun =  bluered ){
	if ( where == 'sample' ){
		n <-as.numeric(as.vector(x$samples[,group] ))
	}else if ( where == 'gene' ) {
		n <-as.numeric(as.vector(x$annotation[,group] ))
	}else {
		stop(paste("Sorry where =",where,"is not supported (only sample and gene)") )
	}
	m <- min( n )
	brks= c( (m-.1),m ,as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
	brks = unique(as.numeric(sprintf("%2.6e", brks)))
	d  <- factor(brks [cut( n, breaks= brks)], levels=brks)
	newN = paste(sep="_",group, 'bins')
	if ( where == 'sample' ){
		x$samples[, newN] <- d
	}else {
		x$annotation[, newN] <- d
	}
	# defined the color
	x$usedObj$colorRange[[newN]] <- colFun( 13 )
	invisible(x)
}