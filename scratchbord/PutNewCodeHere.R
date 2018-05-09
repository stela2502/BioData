addCellCyclePhase <- function( x ) {
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
	
	s.genes <- rownames(x$dat)[which(is.na(match(tolower(rownames(x$dat)),tolower(Seurat::cc.genes$s.genes)))==F)]
	g2m.genes <- rownames(x$dat)[which(is.na(match(tolower(rownames(x$dat)),tolower(Seurat::cc.genes$g2m.genes)))==F)]
	old_m <- ncol(object@meta.data) + 1
	CellCycleScoring(object, g2m.genes, s.genes)
	x$samples <- cbind( x$samples, object@meta.data[, old_m:ncol(object@meta.data)])
	invisible(x)
}
