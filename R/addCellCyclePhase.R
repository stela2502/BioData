#' @name addCellCyclePhase
#' @aliases addCellCyclePhase,BioData-method
#' @rdname addCellCyclePhase-methods
#' @docType methods
#' @description Uses Seurat to address the cell cycle phase
#' @param x the BioData object
#' @param s.genes the s genes (defualt NULL means Seurat::cc.genes$s.genes )
#' @param g2m.genes the s genes (defualt NULL means Seurat::cc.genes$g2m.genes )
#' @title description of function addCellCyclePhase
#' @export 
setGeneric('addCellCyclePhase', ## Name
	function ( x, s.genes = NULL, g2m.genes = NULL ) { ## Argumente der generischen Funktion
		standardGeneric('addCellCyclePhase') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('addCellCyclePhase', signature = c ('BioData'),
	definition = function ( x, s.genes = NULL, g2m.genes = NULL ) {
	if (!requireNamespace("Seurat", quietly = TRUE)) {
		stop("Seurat needed for this function to work. Please install it.",
				call. = FALSE)
	}
	object <- CreateSeuratObject(
			as.matrix(x$dat),
			project = x$name,
			min.cells =3,
			min.genes =1000,
			normalization.method = "LogNormalize",
			scale.factor = 10000
	)
	if ( is.null( s.genes ) ){
		s.genes <- Seurat::cc.genes$s.genes
	}
	if ( is.null( g2m.genes ) ){
		g2m.genes <- Seurat::cc.genes$g2m.genes
	}
	s.genes <- rownames(x$dat)[which(is.na(match(tolower(rownames(x$dat)),tolower(s.genes)))==F)]
	g2m.genes <- rownames(x$dat)[which(is.na(match(tolower(rownames(x$dat)),tolower(g2m.genes)))==F)]
	
	old_m <- ncol(object@meta.data) + 1
	CellCycleScoring(object, g2m.genes, s.genes)
	x$samples <- cbind( x$samples, object@meta.data[, old_m:ncol(object@meta.data)])
	detach( 'Seurat')
	invisible(x)
} )
