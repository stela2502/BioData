#' @name Split4Geo
#' @aliases Split4Geo,BioData-method
#' @rdname Split4Geo-methods
#' @docType methods
#' @description This function will store each column in the data slot table as separate file rbound to the rownames of the data object.
#' The name of the outfile will be the column name + .xls and it will be a tab separated table.
#' An additional md5sums.txt file will be created containing all md5sums for the exported files.
#' @param x the BioData object
#' @param opath the outpath if different from x$outpath
#' @title description of function Split4Geo
#' @export 
if ( ! isGeneric('Split4Geo') ){ setGeneric('Split4Geo', ## Name
	function (x, opath = x$outpath ) { 
		standardGeneric('Split4Geo')
	}
)
}else {
	print ("Onload warn generic function 'Split4Geo' already defined - no overloading here!")
}

setMethod('Split4Geo', signature = c ('BioData'),
	definition = function (x, opath = x$outpath ) {
	if ( ! file.exists( file.path(opath, 'md5sums.txt'))) {
		system(paste( "touch", file.path(opath, 'md5sums.txt') ))
	}
	for ( i in 1:ncol(x$data()) ) {
		d <- cbind( rownames(x$data()), x$data()[,i] )
		colnames(d) <- c( x$rownamescol, colnames(x$dat)[i] )
		write.table( d, file= file.path( opath, paste(colnames(x$dat)[i], 'xls',sep='.' )), sep="\t",quote=F, row.names=F )
		system( paste("md5sum", file.path( opath, paste(colnames(x$dat)[i], 'xls',sep='.' )) , " >>", file.path(opath, 'md5sums.txt') ) )
	}
	invisible(x)
} )
