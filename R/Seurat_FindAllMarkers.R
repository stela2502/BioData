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
	definition = function ( x , condition, logfc.threshold = 1, minPct=0.1 ) {
	
	object <- as_Seurat( x, group=condition, fromRaw=F ) ## I want to use MY data
	
	stats = Seurat::FindAllMarkers( object, logfc.threshold = logfc.threshold, minPct = minPct )

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
	
	x$stats[[name]] = stats
	
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


			stats = cbind( stats, 'p_val_adj' = stats::p.adjust( stats[,'p.value'], method='BH' ) )
			name = paste(sep="_", "Cpp_FindAllMarkers", condition) 
			x$stats[[name]] = stats
			
			invisible(x)
		} )
