#' @name as_Seurat
#' @aliases as_Seurat,BioData-method
#' @rdname as_Seurat-methods
#' @docType methods
#' @description convert BioData to Seurat object but keeping the BioData normalized information.
#' @param x The BioData object
#' @param group Which group should be used as identity in the seurat object (default=NULL)
#' @param fromRaw create the Seurat object from the raw unnormalized data (default =TRUE)
#' @title convert BioData to Seurat
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
		mr = mr[which( ! is.na(mr) ) ]
		mc = match(colnames(x$dat), colnames(object@data))
		mc = mc[which( ! is.na(mc))]
		tmp = x$data()
		bad = which(tmp@x == -1)
		if ( length(bad ) > 0 ){
			tmp@x[bad] = 0
		}
		object@data = tmp[mr,mc]
		rm(tmp)
	}
	if ( ! is.null(group) ) {
		if ( ! is.null(x$samples[,group])) {
			object <- Seurat::SetIdent(object, cells.use = colnames(x$dat), ident.use = as.vector(x$samples[,group]) )
		}
	}
	gc()
	message("Seurat object ready")
	object
} )
