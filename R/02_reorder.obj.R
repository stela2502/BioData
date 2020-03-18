
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
		function ( x, column, ... ) { 
			standardGeneric('reorder.samples')
		}
)
}else {
	print ("Onload warn generic function 'reorder.samples' already defined - no overloading here!")
}

setMethod('reorder.samples', signature = c ('BioData'),
		definition = function ( x, column, ... ) {
			stop("please use the new R6 function instead - more memory efficient!\nobj$reorder.samples(columnName) ")
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
		function ( x, column, ... ) { 
			standardGeneric('reorder.genes')
		}
)
}else {
	print ("Onload warn generic function 'reorder.genes' already defined - no overloading here!")
}
			

setMethod('reorder.genes', signature = c ('BioData'),
		definition = function ( x, column, ... ) {
			stop("please use the new R6 function instead - more memory efficient!\nobj$reorder.genes(columnName) ")
		}

)

