#' @name SeuratVarGenes
#' @aliases SeuratVarGenes,BioData-method
#' @rdname SeuratVarGenes-methods
#' @docType methods
#' @description Pathetic function to use the Seurat FindVariableFeatures function
#' @param x  The BioData object
#' @param n the maount of genes you want (default 300)
#' @title description of function SeuratVarGenes
#' @export 
setGeneric('SeuratVarGenes', ## Name
	function ( x, n=300 ) { ## Argumente der generischen Funktion
		standardGeneric('SeuratVarGenes') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('SeuratVarGenes', signature = c ('BioData'),
	definition = function ( x, n=300  ) {
	m = x$dat
	m@x[which(m@x == -1)] = 0
	var = Seurat::FindVariableFeatures (m)
	x$annotation$VarGeneValue = var[,4]
	x$annotation$VarGene = FALSE
	OK = order(x$annotation$VarGeneValue, decreasing=T)[1:n]
	x$annotation$VarGene [ match( rownames(m)[OK], rownames(x))] = TRUE
	sort(rownames(x)[which(x$annotation$VarGene ==TRUE)])
} )
