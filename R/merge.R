#' @name merge
#' @aliases merge,BioData-method
#' @rdname merge-methods
#' @docType methods
#' @description merges two BioData objects
#' @param x  TEXT MISSING
#' @param objects  TEXT MISSING default=list()
#' @param anno_merge  TEXT MISSING
#' @param anno_uniq  TEXT MISSING
#' @title description of function merge
#' @export 
setGeneric('merge', ## Name
	function ( x, objects=list(), anno_merge, anno_uniq ) { ## Argumente der generischen Funktion
		standardGeneric('merge') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('merge', signature = c ('BioData'),
	definition = function ( x, objects=list(), anno_merge, anno_uniq ) {
	objects = c( x, objects)
	merged <- as_BioData(matrix( 0, ncol=2, nrow=2) )
	gnames = unique(unlist(lapply( 1:length(objects), function( n ) {
								x = objects[[n]]
								if ( is.null(x$samples$sname)){
									x$samples$sname = x$name
								}
								rownames(x$dat)
							} )))
	merged$samples <- objects[[1]]$samples
	lapply(  2:length(objects), function(n) {
				x <- objects[[n]]
				merged$samples <- rbind( merged$samples, x$samples)
			}
	)
	
	dat <- matrix( 0, ,ncol=nrow(merged$samples), nrow=length(gnames))
	rownames(dat) = gnames
	colnames(dat) <- paste(merged$samples$cell_id, merged$samples$sname)
	merged$samples$names = colnames(dat)
	merged$sampleNamesCol <- 'names'
	merged$rownamescol <- A$dat$rownamescol
	
	## and now populate the data frame
	min = 0
	
	for ( n in 1:length(objects) ) {
		
		x <- objects[[n]]
		if ( is.null(x$raw) ) {
			m <- match ( rownames(x$dat), gnames )
			for ( i in 1:ncol(x$dat) ) {
				dat[m,min+i] = x$dat[,i]
			}
			min = min + ncol(x$dat)
		}else {
			m <- match ( rownames(x$raw), gnames )
			for ( i in 1:ncol(x$raw) ) {
				dat[m,min+i] = x$raw[,i]
			}
			min = min + ncol(x$raw)
		}
	}
	
	merged$dat <- as.data.frame(dat)
	
	## and now the annotation table - that is a little bit a problem....
	anno <- matrix( 0, nrow=nrow(dat), ncol=length(anno_merge))
	colnames(anno) <- anno_merge
	
	for (  n in 1:length(objects) ){
		x <- objects[[n]]
		m <-  match( rownames(x$dat), rownames(dat) )
		for ( merge in anno_merge ) {
			anno[m, merge] = as.vector(x$annotation[,merge])
		}
		for ( add in anno_uniq ) {
			new = ncol(anno) +1
			anno <- cbind( anno, rep(0, nrow(anno)) )
			anno[m,new] = as.vector(x$annotation[,add])
			colnames(anno)[new] = paste( sep='_', add, x$name)
		}
		
	}
	merged$annotation = as.data.frame(anno)
	class(merged) <- class(objects[[1]])
	merged
} )
