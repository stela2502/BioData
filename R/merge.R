#' @name merge
#' @aliases merge,BioData-method
#' @rdname merge-methods
#' @docType methods
#' @description merges two or more BioData objects
#' @param x the first BioData object
#' @param objects a list of all others to merge to
#' @title merge 2 or more BioData objects
#' @examples 
#' objects # is a list of BioDat objects having the same gene annotation type (e.g. ENSEMBL IDs!)
#' data <- objects[[1]]
#' objects[[1]] = NULL
#' merged <- merge ( data, objects)
#' @export 
setGeneric('merge', ## Name
	function ( x, objects=list() ) { ## Argumente der generischen Funktion
		standardGeneric('merge') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('merge', signature = c ('BioData'),
	definition = function ( x, objects=list() ) {
		merged = x
	objects = lapply(c(x, objects), function(x){if(is.null(x$rawData)) { x= renew(x) }; x } )
	
	#merge.Matrix( ldat[[paste(n)]], ldat1[[paste(n)]], rownames(ldat[[paste(n)]]), rownames(ldat1[[paste(n)]])  )
	
	
	### checks and fixes in all objects
		
	lapply( objects , function(obj) { 
				merged$dat = Matrix.utils::merge.Matrix( 
						merged$rawData(), 
						obj$rawData(), 
						by.x=rownames(merged$rawData()), 
						by.y=rownames(obj$rawData())
				)
			})
	
	merged$zscored = merged$raw = NULL
	
	lapply( objects , function(obj) { 
				a= t(merged$samples)
				b = t(obj$samples)
				merged$samples = t( Matrix.utils::merge.Matrix( 
						a, 
						b, 
						by.x=rownames(a), 
						by.y=rownames(b)
				))
				colnames(merged$samples) = str_replace_all( colnames(merged$samples), '^y', obj$name )

				a= t(merged$annotation)
				b = t(obj$annotation)
				merged$annotation = t( Matrix.utils::merge.Matrix( 
						a, 
						b, 
						by.x=rownames(a), 
						by.y=rownames(b)
				))
				colnames(merged$annotation) = str_replace_all( colnames(merged$annotation), '^y', obj$name )
			})
	
	class(merged) <- class(objects[[1]])
	gc()
	merged
} )
