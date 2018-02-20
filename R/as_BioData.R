#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @description create a NGSexpressionSet from a counts list object
#' @param dat the counts list you get from Rsubread::featureCounts()
#' @title description of function as_BioData
#' @export 
if ( ! isGeneric('as_BioData') ){ setGeneric('as_BioData', ## Name
	function ( dat, ... ) { ## Argumente der generischen Funktion
		standardGeneric('as_BioData') ## der Aufruf von standardGeneric sorgt für das_BioData Dispatching
	}
)
}else {
	print ("Onload warn generic function 'as_BioData' already defined - no overloading here!")
}

setMethod('as_BioData', signature = c ('list'),
	definition = function ( dat ) {
	ret = NULL
	if (all.equal( names ( dat), c("counts" ,"annotation", "targets", "stat")  ) ) {
		samples <- data.frame(t(dat$stat))
		colnames(samples) <- as.character(t(samples[1,]))
		samples <- samples[-1,]
		samples$filename <- rownames(samples)
		rownames(samples) <-1:nrow(samples)
		ret <- BioData$new( 
				dat= cbind(dat$annotation, dat$counts), 
				Samples = samples, 
				namecol= 'filename', 
				namerow= 'GeneID',
				outpath= ''
		)
	}
	else {
		print ("The list needs to contain the entries counts ,annotation, targets and stat" )
	}
	ret
} )

setMethod('as_BioData', signature = c ('cellexalvr'),
		definition = function ( dat ) {
			dat <- cellexalvr::renew(dat)
			#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
			ok = which(lapply(colnames(dat$meta.cell) , function(x) { all.equal( as.character(as.vector(dat$meta.cell[,x])), colnames(dat$data)) == T } )==T)
			if ( length(ok) == 0) {
				if (all.equal( rownames(dat$meta.cell), colnames(dat$data)) ){
					dat$meta.cell = cbind( cell.name = colnames(dat$data), dat$meta.cell)
					namecol = 'cell.name'
				}
			}
			else {
				namecol = colnames(dat$meta.cell)[ok]
				namecol = namecol[1]
			}
			namerow = NULL
			if (nrow(dat$meta.gene)==0) {
				dat$meta.gene <- matrix(ncol=2, c(rownames(dat$data), rep( 0, nrow(dat$data)) ) )
				colnames(x@meta.gene) = c('Gene.Symbol', 'useless')
				rownames(x@meta.gene) = rownames(x@data)
				namerow = 'Gene.Symbol'
			}else {
				ok = which(lapply(colnames(dat$meta.gene) , function(x) { all.equal( as.character(as.vector(dat$meta.gene[,x])), rownames(dat$data)) == T } )==T)
				if ( length(ok) == 0) {
					if (all.equal( rownames(dat$meta.gene), rownames(dat$data)) ){
						dat$meta.gene = cbind( gene.name = rownames(dat$data), dat$meta.gene)
						namerow = 'gene.name'
					}
				}
				else {
					namerow = colnames(dat$meta.gene)[ok]
					namerow = make.names(namerow[1])
				}
			}
			storage.mode(dat$data) <- 'numeric'
			d <- data.frame(cbind( dat$meta.gene, dat$data))
			if ( nrow(dat$userGroups) == nrow(dat$meta.cell)){
				samples <- data.frame(cbind(dat$meta.cell, dat$userGroups))
			}else {
				samples <- data.frame(dat$meta.cell)
			}
			samples[,namecol] <- make.names(samples[,namecol])
			ret <- BioData$new( d, Samples=samples, name= 'from.cellexalvr', namecol= namecol, namerow=namerow, outpath='./' )
			ret$usedObj <- dat$usedObj
			ret
		}
		
)


setMethod('as_BioData', signature = c ('seurat'),
				definition = function ( dat ) {
					#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
					ok = which(lapply(colnames(dat@meta.data) , function(x) { all.equal( as.character(as.vector(dat@meta.data[,x])), colnames(dat@data)) == T } )==T)
					if ( length(ok) == 0) {
						if (all.equal( rownames(dat@meta.data), colnames(dat@data)) ){
							dat@meta.data = cbind( cell.name = colnames(dat@data), dat@meta.data)
							namecol = 'cell.name'
						}
					}
					else {
						namecol = colnames(dat@meta.data)[ok]
						namecol = namecol[1]
					}
					namerow = NULL
					if (nrow(dat@hvg.info)==0) {
						dat@hvg.info <- matrix(ncol=2, c(rownames(dat@data), rep( 0, nrow(dat@data)) ) )
						colnames(dat@hvg.info) = c('Gene.Symbol', 'useless')
						rownames(dat@hvg.info) = rownames(x@data)
						namerow = 'Gene.Symbol'
					}else {
						ok = which(lapply(colnames(dat@hvg.info) , function(x) { all.equal( as.character(as.vector(dat@hvg.info[,x])), rownames(dat@data)) == T } )==T)
						if ( length(ok) == 0) {
							if (all.equal( rownames(dat@hvg.info), rownames(dat@data)) == TRUE ){
								dat@hvg.info = cbind( gene.name = rownames(dat@data), dat@hvg.info)
								namerow = 'gene.name'
							}else {
								m <- match( rownames(dat@hvg.info), rownames(dat@data))
								dat@hvg.info = cbind( gene.name = rownames(dat@data), dat@hvg.info[order(m),])
								namerow = 'gene.name'
							}
						}
						else {
							namerow = colnames(dat@hvg.info)[ok]
							namerow = make.names(namerow[1])
							## now cast the hvg.info into the same order as the data object.
							m <- match( rownames(dat@hvg.info), rownames(dat@data))
							dat@hvg.info = cbind( gene.name = rownames(dat@data), dat@hvg.info[order(m),])
						}
					}
					#storage.mode(dat$data) <- 'numeric'
					dataa <- as.matrix(dat@data)
					d <- data.frame(cbind( dat@hvg.info[1:2,], dataa[1:2,1:2]))				
					samples <- data.frame(as.matrix(dat@meta.data))
					
					samples[,namecol] <- make.names(samples[,namecol])
					ret <- BioData$new( d, Samples=samples[1:2,], name= 'from.cellexalvr', namecol= namecol, namerow=namerow, outpath='./' )
					
					ret$dat <- dataa
					ret$samples <- samples
					ret$annotation <- dat@hvg.info
					
					ret$zscored <- as.data.frame(as.matrix(dat@scale.data))
					m <- match(colnames(dat@data), colnames(dat@raw.data))
					ret$raw <- as.data.frame(as.matrix(dat@raw.data))[,m]
					## now I need to manually remove unreliable data from the zscored information (set to -20)
					if ( ! is.null( ret$zscored )) {
						d <- function( i, obj) {
							drop <- which(obj$raw[,i] == 0)
							if ( length(drop) > 0 ){
								print ( paste( 'line', i, 'drop', length(drop)) )
								ret$zscored[drop,i] = -20
							}
							0
						}
						lapply(  1:ncol(ret$dat), d, ret)
					}
					
					ret
				}  
)


fetch_first <- function ( sth ) {
	ret <- RSQLite::fetch(sth)
	if (RSQLite::dbHasCompleted(sth)) {
		RSQLite::dbClearResult(sth)
		rm(sth)
	}
	else {
		warn ("query was not finished")
		RSQLite::dbClearResult(sth)
		rm(sth)
	}
	ret
}


SQLite_2_matrix <- function ( fname, useS=NULL, useG=NULL ) {
	dbh <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=fname )
	
	ncol=as.numeric( fetch_first( RSQLite::dbSendQuery(dbh,  "select max(id) from samples" ) ) )
	if ( ! is.null(useS) ) {
		ncol = length(useS)
	}
	nrow=as.numeric( fetch_first( RSQLite::dbSendQuery(dbh,  "select max(id) from genes" ) ) )
	if ( ! is.null(useG) ) {
		ncol = length(useG)
	}
	print ( paste("I create a", ncol,"columns and",nrow,"rows wide matrix"))
	ret <- matrix( 0, nrow=nrow, ncol=ncol )
	
	q = "select gene_id, sample_id, value from datavalues where sample_id = :x"
	
	if ( ! is.null(useG) ) {
		q = paste( q, "and gene_id  IN ( ", paste( collapse=", ", useG),")" )
	}
	sth <- RSQLite::dbSendQuery(dbh, q )
	if(is.null(useS)){
		useS = 1:ncol(ret)
	}
	id = 1
	if (  is.null(useG) ){
		useG = 1:nrow(ret)
	}
	
	
	steps = ceiling(ncol/100)
	pb <- progress_estimated(100)
	print ( "populating the matrix" )
	
	for ( i in  useS ) {
		RSQLite::dbBind(sth, param = list(x = i))
		t <- RSQLite::dbFetch(sth)
		
		m <- match( t$gene_id, useG )
		ret[m, id] <- t$value
		
		if ( id %% steps == 0 ) {
			pb$tick()$print()
			#print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
		}
		id = id +1
	}
	print ( "adding gene and sample names")
	RSQLite::dbClearResult(sth)
	rm(sth)
	q <- "select gname from genes"
	if ( ! is.null(useG) ) {
		q <- paste( q, "where id IN (", paste(collapse=", ", useG),")")
	}
	sth <- RSQLite::dbSendQuery(dbh,q)
	rownames(ret) <- as.character(t(RSQLite::dbFetch(sth)))
	RSQLite::dbClearResult(sth)
	rm(sth)
	
	q <- "select sname from samples"
	if ( ! is.null(useS) ) {
		q <- paste( q, "where id IN (", paste(collapse=", ", useS),")")
	}
	sth <- RSQLite::dbSendQuery(dbh, q )
	colnames(ret) <- as.character(t(RSQLite::dbFetch(sth)))
	RSQLite::dbClearResult(sth)
	RSQLite::dbDisconnect(dbh)
	print ( "done creating the matrix")
	ret
}

setMethod('as_BioData', signature = c ('character'),
		definition = function ( dat, minUMI=100, minGexpr=NULL ) {
			path = dat
			if ( ! file.exists( path ) ) {
				stop( "please - I need an existsing sqlite db at start up!")
			}
			print ( "Obtaining the sample UMI summary")
			sample_UMIs <- SQLite_SampleSummary( dat )
			useS = sample_UMIs$sample_id[which(sample_UMIs$count > minUMI ) ]
			
			useG = NULL
			gene_UMI <- SQLite_ExpressionSummary( dat )
			if ( ! is.null(minGexpr)) {
				gene_UMI <- gene_UMI[which(gene_UMI$'count(value)' > minGexpr )]
				useG = gene_UMI$gene_id
			}
			mat = SQLite_2_matrix ( dat, useS = useS, useG = useG )
			print ( "creating the BioData object" )
			mat <- mat[ -which( apply(mat,1, function(x) { length(which( x != 0)) } ) == 0),]
			#mat$'ensembl_id' <- rownames(mat)
			gene_UMI <- gene_UMI[ match(rownames(mat),gene_UMI$gname ),]
			colnames(mat) <- make.names(colnames(mat))
			ret <- BioData$new( data.frame(cbind(gene_UMI, as.matrix(mat))) ,
					Samples= data.frame( 'cell_id' = colnames(mat) ), namecol='cell_id', namerow= 'gname', name='from.cellranger' )
			class(ret) <- c( 'SingleCells', class(ret))
			ret
		}
)