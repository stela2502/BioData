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
	objects = c( x, objects)
	merged <- as_BioData(matrix( 1, ncol=2, nrow=2, dimnames= list( c('A','B'), c('C','D')) ))
	
	### checks and fixes in all objects
	gnames = unique(unlist(lapply( 1:length(objects), function( n ) {
								x = objects[[n]]
								if ( is.null(x$samples$sname)){
									x$samples$sname = x$name
								}
								if ( is.null( x$samples$nUMI) | is.null(x$samples$nGene) ) {
									if ( ! is.null(x$raw) ){
										x$samples$nUMI <- apply( x$raw,2, sum)
										x$samples$nGene <- apply( x$raw,2, function(a){ length(which(a>0)) })
									}else {
										x$samples$nUMI <- apply( x$dat,2, sum)
										x$samples$nGene <- apply( x$dat,2, function(a){ length(which(a>0)) })
									}
								}
								
								if ( is.null(x$annotation$nUMI) | is.null(x$annotation$nCell) ){
									if ( ! is.null(x$raw) ){
										x$annotation$nUMI <- apply( x$raw,1, sum)
										x$annotation$nCell <- apply( x$raw,1, function(a){ length(which(a>0)) })
									}else {
										x$annotation$nUMI <- apply( x$dat,1, sum)
										x$annotation$nCell <- apply( x$dat,1, function(a){ length(which(a>0)) })
									}
								}
								rownames(x$dat)
							} )))
	
	scols <- colnames(objects[[1]]$samples)
	rcols <- colnames(objects[[1]]$annotation)
	
	snameCol = objects[[1]]$sampleNamesCol
	merged$sampleNamesCol = snameCol
	merged$rownamescol <- x$rownamescol
	
	for ( i in 2:length(objects) ) {
		scols <- intersect( scols, colnames(objects[[i]]$samples))
		rcols <- intersect( rcols, colnames(objects[[i]]$annotation))
		if ( ! snameCol == objects[[i]]$sampleNamesCol){
			stop(paste("the sampleNamesCol of object",(i-1), "named", objects[[i]]$name,"is not the same as in the main object (",snameCol,") - please change that" ) )
		}
		if ( ! x$rownamescol == objects[[i]]$rownamescol){
			stop(paste("the rownamescol of object",(i-1), "named", objects[[i]]$name,"is not the same as in the main object (",x$rownamescol,") - please change that" ) )
		}
	}
	
	merged$samples <- objects[[1]]$samples[,scols]
	lapply(  2:length(objects), function(n) {
				x <- objects[[n]]
				merged$samples <- rbind( merged$samples, x$samples[,scols] )
				NULL
			}
	)
	
	dat <- Matrix( 0 ,ncol=nrow(merged$samples), nrow=length(gnames))
	rownames(dat) = gnames
	colnames(dat) <- make.names(merged$samples[, snameCol])
	merged$samples[, snameCol] <- colnames(dat)
	
	merged$dat <- dat
	
	## The annotation table is not unique for each dataset but should overlap
	## Therefore (1) create the consensus (performed first step in the function. gnames)
	
	merged$annotation <- data.frame( A = gnames)
	colnames(merged$annotation) = x$rownamescol
	
	for( rcolName in setdiff(rcols, x$rownamescol ) ) {
		new_entries = rep(NA, length(gnames) )
		merge=TRUE
		for ( i in 1:length(objects) ) {
			## I need to know if all entries in the datasets do overlap
			m <- match( rownames(objects[[i]]$dat), rownames(merged$dat))
			
			#all.equal(rownames(merged$dat)[m], rownames(objects[[i]]$dat) )
			testable <- intersect(which(!is.na(new_entries)), m)
			if ( length( testable ) > 0 ) {
				t_in_new <- match( rownames(merged$dat)[testable], rownames(objects[[i]]$dat) )
				if ( ! isTRUE( all.equal( 
								as.vector(objects[[i]]$annotation[t_in_new, rcolName]), 
								new_entries[testable]
						) ) 
						){
					merge = FALSE 
				}
			}
			new_entries[ m ] = objects[[i]]$annotation[, rcolName]
		}
		if ( merge ) {
			merged$annotation$'AAADDDEEE' <- new_entries
			colnames(merged$annotation)[ match( 'AAADDDEEE' , colnames(merged$annotation))] = rcolName
		}else {
			### shit - new columns for every entry...
			for ( i in 1:length(objects) ) {
				new_entries <- rep(NA, length(gnames) )
				m <- match( rownames(objects[[i]]$dat), rownames(merged$dat))
				new_entries[m] = objects[[i]]$annotation[, rcolName]
				merged$annotation$'AAADDDEEE' <- new_entries
				colnames(merged$annotation)[ match( 'AAADDDEEE' , colnames(merged$annotation))] = paste(rcolName, objects[[i]]$name )    
			}
		} 
	}	
	
	## and now populate the data frame
	print (paste( "Populate the new dat slot with dim", paste(dim(merged$dat),collapse=":" ) ))
	
	for ( n in 1:length(objects) ) {
		x <- objects[[n]]
		pb <- progress_estimated(100)
		steps = ceiling(ncol(x$dat)/100)
		print ( paste( "Add",ncol(x$dat), "cells from", x$name ))
		if ( is.null(x$raw) ) {
			m <- match ( rownames(x$dat), gnames )
			for ( i in 1:ncol(x$dat) ) {
				mc <- match(make.names(colnames(x$dat)[i]), colnames(merged$dat))
				if ( is.na(mc) ) {
					stop(paste( "Library error - sname" , colnames(x$dat)[i], "not defined in the data colnames", mc ))
				}
				merged$dat[m, mc ] = x$dat[,i]
				if ( i %% steps == 0 ) {
					pb$tick()$print()
					#print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
				}
			}
		}else {
			m <- match ( rownames(x$raw), gnames )
			for ( i in 1:ncol(x$raw) ) {
				mc <- match(make.names(colnames(x$raw)[i]), colnames(merged$dat))
				if ( is.na(mc) ) {
					browser()
					stop(paste( "Library error - sname" , colnames(x$raw)[i], "not defined in the data colnames", mc ))
				}
				merged$dat[m, mc ] = x$raw[,i]
				if ( i %% steps == 0 ) {
					pb$tick()$print()
					#print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
				}
			}
		}
		print ( paste( "added object", x$name ) )
	}
		
	class(merged) <- class(objects[[1]])
	merged
} )
