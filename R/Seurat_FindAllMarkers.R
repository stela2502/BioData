#' @name Seurat_FindAllMarkers
#' @aliases Seurat_FindAllMarkers,BioData-method
#' @rdname Seurat_FindAllMarkers-methods
#' @docType methods
#' @description  Use Seurat's FindAllMarkers function
#' @param x the BioData object
#' @param condition The grouping you wat to analyze
#' @param logfc.threshold test only genes with a liog fold change of default=1.0 or more 
#' @param minPct test only genes with expression in more than default=0.1 fraction of cells in at least one group
#' @title description of function Seurat_FindAllMarkers
#' @export 
setGeneric('Seurat_FindAllMarkers', ## Name
	function ( x , condition, logfc.threshold=1.0, minPct=0.1) { 
		standardGeneric('Seurat_FindAllMarkers')
	}
)

setMethod('Seurat_FindAllMarkers', signature = c ('BioData'),
	definition = function ( x , condition, logfc.threshold=1.0, minPct=0.1) {
	
	
	#object <- as_Seurat( x, group=condition, fromRaw=F ) ## I want to use MY data

	
			
	CppStats <- function( n ) {
		print ( paste("Calc wilcox test for",n) )
		OK = which(grp.vec == n )
		BAD= which(grp.vec != n )
		if ( length(OK) > 0 && length(BAD) > 0 ){
			r = as.data.frame(
					FastWilcoxTest::StatTest( Matrix::t( x$dat), OK, BAD, logfc.threshold, minPct )
			)
			r= r[order( r[,'p.value']),]
			r = cbind( r, cluster= rep(n,nrow(r) ), gene=rownames(x$dat)[r[,1]] )
		}
		r
	}
	
	stats = NULL;
	grp.vec = as.vector(x$samples[,condition])
	
	for ( n in  unique( sort(grp.vec)) ) {
		stats = rbind( stats, CppStats(n) )
	}
	
	
	add_to_stat <- function( x, stat, name ) {
		if ( ! is.na( match( name, names(x$stats)))){
			x$stats[[ match( name, names(x$stats)) ]] <- stat
		}else {
			x$stats[[ length( x$stats ) +1 ]] <- stat
			names(x$stats)[length(x$stats) ] <- name
		}
		x
	}
	
	x <- add_to_stat ( x, stats, paste(sep="_", "wilcox_FindAllMarkers", condition) )
	
	invisible(x)
} )
