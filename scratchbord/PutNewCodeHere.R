convert_to <- function(x, type=c("MAST", "DESeq2", "scran") , species=NULL, ...) {
	ret <- NULL
	toM <- function (x) {
		d <- as.matrix(x)
		d[which(d==-20)] <- NA
		d[is.na(d)] <- 0
		d
	}
	if (type == 'DESeq2' ) {
		ret <- DESeq2::DESeqDataSetFromMatrix(
				toM(x$raw),
				x$samples
		)
	}else if ( type== "MAST" ) {
		if ( is.null(species) ) {
			stop ( "to convert to scran I need a species or AnnotationDbi object" )
			## so now I need to get the ensembl ids into this crappy object
			ret <- scran::SingleCellExperiment(list(counts=x$raw))
			colData(ret) <- x$samples
			rowData(ret) <- x$annotation
		}
	}
	ret
}