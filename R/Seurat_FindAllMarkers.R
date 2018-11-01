#' @name Seurat_FindAllMarkers
#' @aliases Seurat_FindAllMarkers,BioData-method
#' @rdname Seurat_FindAllMarkers-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param condition  TEXT MISSING
#' @title description of function Seurat_FindAllMarkers
#' @export 
setGeneric('Seurat_FindAllMarkers', ## Name
	function ( x , condition) { ## Argumente der generischen Funktion
		standardGeneric('Seurat_FindAllMarkers') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('Seurat_FindAllMarkers', signature = c ('BioData'),
	definition = function ( x , condition) {
	
	
	object <- as_Seurat( x, group=condition, fromRaw=F ) ## I want to use MY data

	stats = Seurat::FindAllMarkers( object )
	
	x <- add_to_stat ( x, stats, paste(sep="_", "Seurat_FindAllMarkers", condition) )
	
	invisible(x)
} )
