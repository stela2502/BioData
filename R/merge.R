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
	function ( x, objects=list() ) { 
		standardGeneric('merge')
	}
)

setMethod('merge', signature = c ('BioData'),
		definition = function ( x, objects=list() ) {
			objects = lapply(c(x, objects), function(x){if(is.null(x$rawData)) { x= renew(x) }; x } )
			
			merged <- as_BioData(matrix( 1, ncol=2, nrow=2, dimnames= list( c('A','B'), c('C','D')) ))
			
			#merge.Matrix( ldat[[paste(n)]], ldat1[[paste(n)]], rownames(ldat[[paste(n)]]), rownames(ldat1[[paste(n)]])  )
			
			
			### checks and fixes in all objects
			gnames = unique(unlist(lapply( 1:length(objects), function( n ) {
										x = objects[[n]]
										if ( is.null(x$samples$sname)){
											x$samples$sname = x$name
										}
										if ( is.null( x$samples$nUMI) | is.null(x$samples$nGene) ) {
											if ( ! is.null(x$raw) ){
												x$samples$nUMI <- Matrix::colSums( x$raw )
												x$samples$nGene <- apply( x$raw,2, function(a){ length(which(a>0)) })
											}else {
												x$samples$nUMI <- Matrix::colSums( x$dat )
												x$samples$nGene <- apply( x$dat,2, function(a){ length(which(a>0)) })
											}
										}
										
										if ( is.null(x$annotation$nUMI) | is.null(x$annotation$nCell) ){
											if ( ! is.null(x$raw) ){
												x$annotation$nUMI <- Matrix::rowSums( x$raw )
												x$annotation$nCell <- apply( x$raw,1, function(a){ length(which(a>0)) })
											}else {
												x$annotation$nUMI <- Matrix::rowSums( x$dat )
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
			
			dat <- Matrix::Matrix( 0 ,ncol=nrow(merged$samples), nrow=length(gnames))
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
			
			values= NULL;

			## adding to Matrix - slow as hell
			## pre_defining the vectors (adding to them) - slow as hell
			## creating a temporary file - a huge improvement, but a combo of the last 2 should be best
			## creating a table if data is more than 10000 ?
			#file = filePath( Sys.getenv('SNIC_TMP'), 'tmp.table.sqlite')
			#if ( file.exists(file)) {unlink( file )}
			#dbh <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=file )
			#dbh <- RSQLite::dbConnect(RSQLite::SQLite(),":memory:" )
			# OK as I anyhow put all this info into memory - why not allocate the vectors according to the requirements?
			
			for ( n in 1:length(objects) ) {
				x <- objects[[n]]
				pb <- dplyr::progress_estimated(100)
				steps = ceiling(ncol(x$dat)/100)
				print ( paste( "Add",ncol(x$dat), "cells from", x$name ))
				
				melted = FastWilcoxTest::meltSparseMatrix( x$rawData() )
				## now I need to translate the ids for genes and cells
				melted[,1] = match( rownames(x$dat), rownames(merged$dat))[melted[,1]]
				melted[,2] = match( colnames(x$dat), colnames(merged$dat))[melted[,2]]
				values = rbind( values, melted )
				
				print ( paste( "added object", x$name ) )
			}
			#browser()
			tmp <- Matrix::sparseMatrix( i=values[,1], j=values[,2], x=values[,3], dims=dim(merged$dat) )
			
			#rm(dbh) 
			colnames(tmp) <- colnames(merged$dat)
			rownames(tmp) <- rownames(merged$dat)
			merged$dat = tmp
			
			rm( tmp, values)
			class(merged) <- class(objects[[1]])
			gc()
			merged
		} )
