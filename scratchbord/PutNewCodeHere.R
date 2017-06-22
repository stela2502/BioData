Split4Geo <- function(x, opath = x$outpath ) {
	if ( ! file.exists( file.path(opath, 'md5sums.txt'))) {
		system(paste( "touch", file.path(opath, 'md5sums.txt') ))
	}
	for ( i in 1:ncol(x$data) ) {
		d <- cbind( rownames(x$data), x$data[,i] )
		colnames(d) <- c( x$rownamescol, colnames(x$data)[i] )
		write.table( d, file= file.path( opath, paste(colnames(x$data)[i], 'xls',sep='.' )), sep="\t",quote=F, row.names=F )
		system( paste("md5sum", file.path( opath, paste(colnames(x$data)[i], 'xls',sep='.' )) , " >>", file.path(opath, 'md5sums.txt') ) )
	}
	invisible(x)
}

