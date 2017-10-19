#' @name defineGOIs
#' @aliases defineGOIs,cellexalvr-method
#' @rdname defineGOIs-methods
#' @docType methods
#' @description  Allows the user to define (G)enes (O)f (I)nterest lists in the object
#' @param x A BioData object
#' @param name the name of the GIO list (eg TFs or epigenetic)
#' @param genes a list of gene symbols that match to the $data rownames
#' @param lables a list of lables for the GIO column (default NULL)
#' @param gene_col which column in the annotation data contains the IDs you give me
#' default == rownames(x$data)
#' @title description of function defineGOIs
#' @export defineGOIs
setGeneric('defineGOIs', ## Name
	function ( x,name, genes, lables=NULL, ...) { ## Argumente der generischen Funktion
		standardGeneric('defineGOIs') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('defineGOIs', signature = c ('BioData'),
	definition = function ( x, name, genes, lables=NULL, gene_col=NULL, ... ) {
		if ( is.null(lables) ) {
			lables = rep(name, length(genes))
		}
		if ( ! is.null(gene_col)) {
			genes = rownames(x$data)[ match(genes, x$annotation[,gene_col]) ]
		}
		if (nrow(x$annotation)==0) {
			x$annotation <- matrix(ncol=2, c(rownames(x$data), rep( 0, nrow(x$data)) ) )
			colnames(x$annotation) = c('Gene Symbol', 'useless')
			rownames(x$annotation) = rownames(x$data)
		}
		if ( ! is.na( match(name, colnames(x$annotation)))) {
			stop( "Sorry, but this GIO list has already been defined" )
		}
		
		#OK <- which(is.na(match( rownames(x$data), genes)) == F)
		OK_genes <- match( genes, rownames(x$data) )
		OK <- OK_genes[which(is.na(OK_genes) ==F)]
		n = rep(NA, nrow(x$annotation))
		#browser()
		n[OK] = as.vector(lables[which(is.na(OK_genes) ==F)])
		x$annotation <- cbind( x$annotation, n )
		colnames(x$annotation)[ncol(x$annotation)] = name
		x$usedObj$GOIs = c( x$usedObj$GOIs, name)
		
		invisible(x)
	} 
)
