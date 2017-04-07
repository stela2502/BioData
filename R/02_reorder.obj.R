
#' @name reorder.samples
#' @aliases reorder.samples,BioData-method
#' @rdname reorder.samples-methods
#' @docType methods
#' @description this function reorderes the BioData object based on a column in the annotation table (e.g. for plotting)
#' @param x the BioData object
#' @param column the annotation column to reorder on
#' @title description of function remove.genes
#' @export 
setGeneric('reorder.samples', ## Name
		function ( x, column, ... ) { ## Argumente der generischen Funktion
			standardGeneric('reorder.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('reorder.samples', signature = c ('BioData'),
		definition = function ( x, column, ... ) {
			x$data <- x$data[ , order( x$samples[,column])]
			x$samples <- x$samples[order( x$samples[,column]),]
			invisible(x)
		}
)

#' @name reorder.genes
#' @aliases reorder.genes,BioData-method
#' @rdname reorder.genes-methods
#' @docType methods
#' @description this function reorderes the BioData object based on a column in the samples table (e.g. for plotting)
#' @param x the BioData object
#' @param column the samples column to reorder on
#' @title description of function remove.genes
#' @export 
setGeneric('reorder.genes', ## Name
		function ( x, column, ... ) { ## Argumente der generischen Funktion
			standardGeneric('reorder.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('reorder.genes', signature = c ('BioData'),
		definition = function ( x, column, ... ) {
			x$data <- x$data[ order( x$annotation[,column]),]
			x$annotation <- x$annotation[order( x$annotation[,column]),]
			invisible(x)
		}

)

