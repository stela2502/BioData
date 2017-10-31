
#' @name reorder.samples
#' @aliases reorder.samples,BioData-method
#' @rdname reorder.samples-methods
#' @docType methods
#' @description this function reorderes the BioData object based on a column in the annotation table (e.g. for plotting)
#' @param x the BioData object
#' @param column the annotation column to reorder on
#' @title description of function remove.genes
#' @export 
if ( ! isGeneric('reorder.samples') ){ 
	setGeneric('reorder.samples', ## Name
		function ( x, column, ... ) { ## Argumente der generischen Funktion
			standardGeneric('reorder.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
}else {
	print ("Onload warn generic function 'reorder.samples' already defined - no overloading here!")
}

setMethod('reorder.samples', signature = c ('BioData'),
		definition = function ( x, column, ... ) {
			idx <-  order( x$samples[,column])
			x$dat <- x$dat[ , idx ]
			if( !is.null(x$zscored)) {
				x$zscored <- x$zscored[ , idx]
			}
			if( !is.null(x$raw)) {
				x$raw <- x$raw[ , idx]
			}
			x$samples <- x$samples[idx ,]
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
if ( ! isGeneric('reorder.genes') ){ setGeneric('reorder.genes', ## Name
		function ( x, column, ... ) { ## Argumente der generischen Funktion
			standardGeneric('reorder.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
}else {
	print ("Onload warn generic function 'reorder.genes' already defined - no overloading here!")
}
			

setMethod('reorder.genes', signature = c ('BioData'),
		definition = function ( x, column, ... ) {
			idx <-  order( x$annotation[,column])
			x$dat <- x$dat[ idx ,]
			if( !is.null(x$zscored)) {
				x$zscored <- x$zscored[ idx, ]
			}
			if( !is.null(x$raw)) {
				x$raw <- x$raw[ idx, ]
			}
			t <- x$annotation[ idx,]
			if ( class(t) == 'factor' ) {
				t <- data.frame(t)
				colnames(t) <- colnames(x$annotation)
			}
			x$annotation <- t
			invisible(x)
		}

)

