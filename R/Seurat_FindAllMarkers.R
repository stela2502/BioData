#' @name Seurat_FindAllMarkers
#' @aliases Seurat_FindAllMarkers,BioData-method
#' @rdname Seurat_FindAllMarkers-methods
#' @docType methods
#' @description  Use Seurat's FindAllMarkers function
#' @param x the BioData object
#' @param condition The grouping you wat to analyze
#' @param logfc.threshold  Seurat FindAllmarkers option used for speed up (default = 1) 
#' @title description of function Seurat_FindAllMarkers
#' @export 
setGeneric('Seurat_FindAllMarkers', ## Name
	function ( x , condition, logfc.threshold = 1 ) { 
		standardGeneric('Seurat_FindAllMarkers')
	}
)

setMethod('Seurat_FindAllMarkers', signature = c ('BioData'),
	definition = function ( x , condition, logfc.threshold = 1 ) {
	
	object <- as_Seurat( x, group=condition, fromRaw=F ) ## I want to use MY data
	
	stats = Seurat::FindAllMarkers( object, logfc.threshold = logfc.threshold )
	
	add_to_stat <- function( x, stat, name ) {
		if ( ! is.na( match( name, names(x$stats)))){
			x$stats[[ match( name, names(x$stats)) ]] <- stat
		}else {
			x$stats[[ length( x$stats ) +1 ]] <- stat
			names(x$stats)[length(x$stats) ] <- name
		}
		x
	}
	
	x <- add_to_stat ( x, stats, paste(sep="_", "Seurat_FindAllMarkers", condition) )
	
	invisible(x)
} )

#' @name Cpp_FindAllMarkers
#' @aliases Cpp_FindAllMarkers,BioData-method
#' @rdname Cpp_FindAllMarkers-methods
#' @docType methods
#' @description  Use Seurat's FindAllMarkers function
#' @param x the BioData object
#' @param condition The grouping you wat to analyze
#' @param logfc.threshold seurat option used for speed up (default = 1)
#' @param  minPct only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' @param onlyPos check only higher expression (default FALSE)	
#' @title description of function Cpp_FindAllMarkers
#' @export 
setGeneric('Cpp_FindAllMarkers', ## Name
		function ( x , condition, logfc.threshold = 1, minPct = 0.1, onlyPos = FALSE) { 
			standardGeneric('Cpp_FindAllMarkers')
		}
)

setMethod('Cpp_FindAllMarkers', signature = c ('BioData'),
		definition = function ( x , condition, logfc.threshold = 1, minPct = 0.1 , onlyPos = FALSE) {
			
			grp.vec = as.vector( x$samples[,condition])
			CppStats <- function( n ) {
				OK = which(grp.vec == n )
				BAD= which(grp.vec != n )
				r = as.data.frame(FastWilcoxTest::StatTest( mat, OK, BAD, logfc.threshold, minPct, onlyPos ))
				r= r[order( r[,'p.value']),]
				r = cbind( r, cluster= rep(n,nrow(r) ), gene=rownames(x$dat)[r[,1]] )
				r
			}
			
			mat = Matrix::t(x$dat)
			bad = which( mat@x <= 0 )
			if ( length(bad) > 0 ){
				mat@x[bad] = 0;
			}
			
			all_markers = NULL;
			for ( n in  unique( sort(grp.vec)) ) {
				all_markers = rbind( all_markers, CppStats(n) )
			}
			all_markers = cbind( all_markers, 'p_val_adj' = stats::p.adjust( all_markers[,'p.value'], method='BH' ) )
			name = paste(sep="_", "Cpp_FindAllMarkers", condition) 
			x$stats[[name]] = all_markers
			
			invisible(x)
		} )