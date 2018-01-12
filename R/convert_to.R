#' @name convert_to
#' @aliases convert_to,BioData-method
#' @rdname convert_to-methods
#' @docType methods
#' @description convert a biodata object into another object type
#' @param x the BioData object
#' @param type the type to convert to 'scran' should be the only working here....
#' @param species if converting to scran you need to tell me which species the data comes from in order to add the ensembol ids scan needs
#' @param ...  TEXT MISSING
#' @title description of function convert_to
#' @export 
setGeneric('convert_to', ## Name
	function (x, type=c("MAST", "DESeq2", "scran") , species=NULL, ...) { ## Argumente der generischen Funktion
		standardGeneric('convert_to') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('convert_to', signature = c ('BioData'),
	definition = function (x, type=c("MAST", "DESeq2", "scran") , species=NULL, ...) {
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
	}
	if ( type== "scran" ) {
		if (!requireNamespace("scran", quietly = TRUE)) {
			stop("scran needed for this function to work. Please install it.",
					call. = FALSE)
		}
		if ( is.null(species) ) {
			stop ( "to convert to scran I need a species or AnnotationDbi object" )
		}
			## so now I need to get the ensembl ids into this crappy object
		if ( is.na( match( 'ensembl_id', colnames(x$annotation))) ) {
			if ( is.na( match( 'gene.name', colnames(x$annotation)))) {
				stop( "I need either a ensembl_id or a 'gene.name' annotation column in order to convert to scran")
			}
			x$annotation$ensembl_id = getGeneInfo( as.vector(x$annotation$gene.name), 
					species = species, from = "symbol", what='ensembl_id', what_tab="ensembl")
		}
		x <- x$clone()		
		changeNames(x, 'row', 'ensembl_id')
		rownames(x$annotation) <- rownames(x$dat)
		rownames(x$samples) <- colnames(x$dat)
		ret <- SingleCellExperiment(list(counts=as.matrix(x$raw)))
		rowData(ret) <- as.matrix(x$annotation)
		colData(ret) <- DataFrame(x$samples)
		if ( ! is.null(x$annotation$gene.name) ) {
			rowData(ret)$gene.symbol <- rowData(ret)$gene.name
			rowData(ret)$gene.name <- rowData(ret)$ensembl_id.unique
		}
	}
	ret
} )
