
Seurat_FindAllMarkers <- function( x , condition){
	
	
	object <- as_Seurat( x, group=condition )

	stats = FindAllMarkers( object )
	
	x <- add_to_stat ( x, stats, paste(sep="_", "Seurat_FindAllMarkers", condition) )
	
	invisible(x)
}