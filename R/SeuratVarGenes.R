#' @name SeuratVarGenes
#' @aliases SeuratVarGenes,BioData-method
#' @rdname SeuratVarGenes-methods
#' @docType methods
#' @description Pathetic function to use the Seurat FindVariableFeatures function
#' @param x  The BioData object
#' @param minVar the min scaled varianze (default=1)
#' @title description of function SeuratVarGenes
#' @export 
setGeneric('SeuratVarGenes', ## Name
	function ( x, minVar=1 ) { ## Argumente der generischen Funktion
		standardGeneric('SeuratVarGenes') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('SeuratVarGenes', signature = c ('BioData'),
	definition = function ( x, minVar=1 ) {
	m = x$dat
	m@x[which(m@x == -1)] = 0
	var = Seurat::FindVariableFeatures (m)
	x$annotation$VarGeneValue = var[,4]
	x$annotation$VarGene = FALSE
	x$annotation$VarGene [ match( rownames(m)[which(var[,4] > minVar)], rownames(x))] = TRUE
	invisible(x)
} )
