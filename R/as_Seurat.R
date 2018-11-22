#' @name as_Seurat
#' @aliases as_Seurat,BioData-method
#' @rdname as_Seurat-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param group  TEXT MISSING default=NULL
#' @title description of function as_Seurat
#' @export 
setGeneric('as_Seurat', ## Name
	function ( x , group=NULL, fromRaw = T ) { 
		standardGeneric('as_Seurat')
	}
)

setMethod('as_Seurat', signature = c ('BioData'),
	definition = function ( x , group=NULL, fromRaw = T ) {
	
	object <- NULL
	message('Create Seurat object from raw data')
	if (  ! is.null(x$raw)) {
		object <- Seurat::CreateSeuratObject(
				x$raw,
				project = x$name,
				min.cells =1,
				min.genes =1,
				normalization.method = "LogNormalize",
				scale.factor = 10000
		)
	}else {
		object <- Seurat::CreateSeuratObject(
				x$dat,
				project = x$name,
				min.cells =1,
				min.genes =1,
				normalization.method = "LogNormalize",
				scale.factor = 10000
		)
	}
	if ( ! fromRaw ) {
		message("replacing Seurat norm data with BioData norm data")
		mr = match(rownames(x$dat), rownames(object@data))
		mr = mr[which( ! is.na(rm) ) ]
		mc = match(colnames(x$dat), colnames(object@data))
		object@data = x$data()[mr,mc]
		bad = which(object@data@x == -1)
		if ( length(bad ) > 0 ){
			object@data@x[bad] =0;
		}
	}
	if ( ! is.null(group) ) {
		if ( ! is.null(x$samples[,group])) {
			object <- Seurat::SetIdent(object, cells.use = colnames(x$dat), ident.use = as.vector(x$samples[,group]) )
		}
	}
	
	object
} )
