#' @name export.data
#' @aliases export.data,BioData-method
#' @rdname export.data-methods
#' @docType methods
#' @description  write the BioData data, annotation and samples tables to disk
#' @param x the BioData object
#' @title description of function write.data
#' @export 
if ( ! isGeneric('export.data') ){ setGeneric('export.data', ## Name
		function ( x ) { 
			standardGeneric('export.data')
		}
)
}else {
	print ("Onload warn generic function 'export.data' already defined - no overloading here!")
}

setMethod('export.data', signature = c ( 'BioData') ,
		definition = function ( x ) {
			ofile = paste(x$outpath,str_replace_all(x$name, "\\s+", "_" ),sep='')
			write.table(cbind ( x$annotation, x$dat ), file=paste(ofile,"expressionValues.xls",sep='_'),sep='\t', row.names=F,quote=F )
			if ( ! is.null(x$raw) ) {
				write.table(cbind ( x$annotation, x$raw ), file=paste(ofile,"raw_expressionValues.xls",sep='_'),sep='\t', row.names=F,quote=F )
			}
			if ( ! is.null( x$zscored) ) {
				write.table(cbind ( x$annotation, x$zscored ), file=paste(ofile,"zscored_expressionValues.xls",sep='_'),sep='\t', row.names=F,quote=F )
			}
			write.table( x$samples, file=paste( ofile, "sampleInformation.xls",sep='_' ), sep='\t', row.names=F,quote=F )
	#		write.table( x$annotation, file=paste( ofile, "gene_annotation.xls",sep='_' ), sep='\t', row.names=F,quote=F )
		})

