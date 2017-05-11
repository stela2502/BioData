#' @name export.data
#' @aliases export.data,BioData-method
#' @rdname export.data-methods
#' @docType methods
#' @description  write the BioData data, annotation and samples tables to disk
#' @param x the BioData object
#' @title description of function write.data
#' @export 
setGeneric('export.data', ## Name
		function ( x ) { ## Argumente der generischen Funktion
			standardGeneric('export.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('export.data', signature = c ( 'BioData') ,
		definition = function ( x ) {
			ofile = paste(x$outpath,str_replace_all(x$name, "\\s+", "_" ),sep='')
			write.table(cbind ( x$annotation, x$data ), file=paste(ofile,"expressionValues.xls",sep='_'),sep='\t', row.names=F,quote=F )
			write.table( x$samples, file=paste( ofile, "sampleInformation.xls",sep='_' ), sep='\t', row.names=F,quote=F )
	#		write.table( x$annotation, file=paste( ofile, "gene_annotation.xls",sep='_' ), sep='\t', row.names=F,quote=F )
		})

