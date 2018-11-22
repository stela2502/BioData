#' @name addCellCyclePhase
#' @aliases addCellCyclePhase,BioData-method
#' @rdname addCellCyclePhase-methods
#' @docType methods
#' @description Uses Seurat to address the cell cycle phase
#' @param x the BioData object
#' @param s.genes the s genes (defualt NULL means Seurat::cc.genes$s.genes )
#' @param g2m.genes the s genes (defualt NULL means Seurat::cc.genes$g2m.genes )
#' @param gnameCol the name of the Gene Symbol column in the annotation table (default NULL) 
#' @title description of function addCellCyclePhase
#' @export 
setGeneric('addCellCyclePhase', ## Name
	function ( x, s.genes = NULL, g2m.genes = NULL, gnameCol=NULL ) { 
		standardGeneric('addCellCyclePhase')
	}
)

setMethod('addCellCyclePhase', signature = c ('BioData'),
	definition = function ( x, s.genes = NULL, g2m.genes = NULL, gnameCol=NULL ) {
	if (!requireNamespace("Seurat", quietly = TRUE)) {
		stop("Seurat needed for this function to work. Please install it.",
				call. = FALSE)
	}
	object <- as_Seurat( x )
	
	if ( is.null( s.genes ) ){
		s.genes <- Seurat::cc.genes$s.genes
	}
	if ( is.null( g2m.genes ) ){
		g2m.genes <- Seurat::cc.genes$g2m.genes
	}
	if ( is.null(gnameCol) ) {
		s.genes <- rownames(x$dat)[which(is.na(match(tolower(rownames(x$dat)),tolower(s.genes)))==F)]
		g2m.genes <- rownames(x$dat)[which(is.na(match(tolower(rownames(x$dat)),tolower(g2m.genes)))==F)]
	}else{
		s.genes <- rownames(x$dat)[which(is.na(match(tolower(as.vector(x$annotation[,gnameCol])),tolower(s.genes)))==F)]
		g2m.genes <- rownames(x$dat)[which(is.na(match(tolower(as.vector(x$annotation[,gnameCol])),tolower(g2m.genes)))==F)]
	}
	if (!( length(s.genes) >0 &  length(s.genes) >0 ) ) {
		stop(paste( "Sorry - I could not find enough S and G2M specific genes in the dataset (need to provide a gnameCol?)" ) )
	}
	old_m <- ncol(object@meta.data) + 1
	message ( "addressing the cell cycle phase using Seurat::CellCycleScoring")
	object <- Seurat::CellCycleScoring(object, g2m.genes, s.genes)
	x$samples <- cbind( x$samples, object@meta.data[, old_m:ncol(object@meta.data)])
	rm(object)
	gc()
	#eval(detach( 'package:Seurat'))
	invisible(x)
} )
