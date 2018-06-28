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
	pI = sum(unlist(lapply(objects, function(x) {
								if ( is.null(x$raw)) {
									length(which(x$dat > 0))
								}else {
									length(which(x$raw > 0))
								}
							})))
	
	pM = pI * 2 
	#pI = 20000 
	r = rep( NA, pI ) #rows
	j = rep( NA, pI) #cols
	v = rep( NA, pI) #vals
	pos = 0
	## adding to Matrix - slow as hell
	## pre_defining the vartors (adding to them) - slow as hell
	## creating a temporary file - a huge imrovement, but a combo of the last 2 should be best
	## creating a table if data is more than 10000 ?
	#file = filePath( Sys.getenv('SNIC_TMP'), 'tmp.table.sqlite')
	#if ( file.exists(file)) {unlink( file )}
	#dbh <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=file )
	#dbh <- RSQLite::dbConnect(RSQLite::SQLite(),":memory:" )
	# OK as I anyhow put all this info into memory - why not allocate the vectors according to the requirements?
    
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
				ok = which(x$dat[,i] > 0)
				if ( length(ok) > 0 ) {
					#write.table( cbind( r= m[ok] , j= rep(mc, length(ok), v= x$dat[ok,i]) ),file=file, quote=F, append=T, row.names=F, col.names=F )
				    add = (pos+1):(pos+length(ok)+1)
					r[add] = m[ok] ## rows
					j[add] = rep(mc, length(ok) ) ## cols
					v[add] = x$dat[ok,i] ## vals
					pos = pos+length(ok)
				}
				if ( i %% steps == 0 ) {
					pb$tick()$print()
					#print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
				}
			}
			if ( pos > pM ) {
				RSQLite::dbWriteTable(dbh,'datavalues',data.frame( r= r[1:pos] , j= j[1:pos], v= v[1:pos] ),append = TRUE)
				#write.table( cbind( r= r[1:pos] , j= j[1:pos], v= v[1:pos] ),file=file, quote=F, append=T, row.names=F, col.names=F )
				r = rep( NA, pI ) #rows
				j = rep( NA, pI) #cols
				v = rep( NA, pI) #vals
				pos = 0
			}
		}else {
			m <- match ( rownames(x$raw), gnames )
			for ( i in 1:ncol(x$raw) ) {
				mc <- match(make.names(colnames(x$raw)[i]), colnames(merged$dat))
				if ( is.na(mc) ) {
					browser()
					stop(paste( "Library error - sname" , colnames(x$raw)[i], "not defined in the data colnames", mc ))
				}
				ok = which(x$raw[,i] > 0)
				if ( length(ok) > 0 ) {
					#write.table( cbind( r= m[ok] , j= rep(mc, length(ok), v= x$raw[ok,i]) ),file=file, quote=F, append=T, row.names=F, col.names=F )
					add = (pos+1):(pos+length(ok)+1)
					r[add] = m[ok] ## rows
					j[add] = rep(mc, length(ok) ) ## cols
					v[add] = x$dat[ok,i] ## vals
					pos = pos+length(ok)
				}
				if ( i %% steps == 0 ) {
					pb$tick()$print()
					#print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
				}
			}
			if ( pos > pM ) {
				RSQLite::dbWriteTable(dbh,'datavalues', data.frame( r= r[1:pos] , j= j[1:pos], v= v[1:pos] ),append = TRUE)
				#write.table( cbind( r= r[1:pos] , j= j[1:pos], v= v[1:pos] ),file=file, quote=F, append=T, row.names=F, col.names=F )
				r = rep( NA, pI ) #rows
				j = rep( NA, pI) #cols
				v = rep( NA, pI) #vals
				pos = 0
			}
		}
		print ( paste( "added object", x$name ) )
	}
	# not required - all written to a vector anyhow!
#	if ( pos > 0 ){ ## clean up
#		RSQLite::dbWriteTable(dbh,'datavalues',data.frame( r= r[1:pos] , j= j[1:pos], v= v[1:pos] ),append = TRUE)
#		#write.table( cbind( r= r[1:pos] , j= j[1:pos], v= v[1:pos] ),file=file, quote=F, append=T, row.names=F, col.names=F )
#	}
	#browser()
	#print ("loading tmp database")
	#q = "select * from datavalues"
	#sth <- RSQLite::dbSendQuery(dbh, q )
	#t <- RSQLite::dbFetch(sth)
	
	tmp <- sparseMatrix( i=r, j=j, x=v, dims=dim(merged$dat) )
	
	#rm(dbh) 
	colnames(tmp) <- colnames(merged$dat)
	rownames(tmp) <- rownames(merged$dat)
	merged$dat = tmp
	
	class(merged) <- class(objects[[1]])
	gc()
	merged
} )
