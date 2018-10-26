#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @description create a NGSexpressionSet from an other object
#' @export 
if ( ! isGeneric('as_BioData') ){ setGeneric('as_BioData', ## Name
	function ( dat, ... ) { ## Argumente der generischen Funktion
		standardGeneric('as_BioData') ## der Aufruf von standardGeneric sorgt für das_BioData Dispatching
	}
)
}else {
	print ("Onload warn generic function 'as_BioData' already defined - no overloading here!")
}

#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @param dat a list coming for Rsubread scanning of NGS data
#' @title description of function as_BioData
#' @export 
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
				dat= Matrix(dat$counts), 
				Samples = samples, 
				annotation= dat$annotation,
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

#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @param dat a matrix object with cells in columns and genes in rows
#' @title description of function as_BioData
#' @export 
setMethod('as_BioData', signature = c ('matrix'),
		definition = function ( dat ) {
				dat <- Matrix(dat)
				snames <- colnames(dat)
				
			ret <- BioData$new( 
					dat=  dat, 
					Samples = data.frame('sampleName' = snames ), 
					annotation= data.frame( GeneID = rownames(dat) ),
					namecol= 'sampleName', 
					namerow= 'GeneID',
					outpath= ''
			)
			ret
		}
)

setMethod('as_BioData', signature = c ('data.frame'),
		definition = function ( dat ) {
			snames <- colnames(dat)
			#dat$GeneID = rownames(dat)
			ret <- BioData$new( 
					dat=  Matrix(as.matrix(dat)), 
					Samples = data.frame('sampleName' = snames ), 
					annotation= data.frame( 'GeneID' = rownames(dat)),
					namecol= 'sampleName', 
					namerow= 'GeneID',
					outpath= ''
			)
			ret
		}
)

#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @param dat a cellexlvrR object
#' @title description of function as_BioData
#' @export 
setMethod('as_BioData', signature = c ('cellexalvr'),
		definition = function ( dat ) {
			#dat <- cellexalvr::renew(dat)
			#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
			ok = which(lapply(colnames(dat@meta.cell) , function(x) { all.equal( as.character(as.vector(dat@meta.cell[,x])), colnames(dat@data)) == T } )==T)
			namecol = NULL
			if ( length(ok) == 0) {
				if (all.equal( rownames(dat@meta.cell), colnames(dat@data)) == TRUE ){
					dat@meta.cell = cbind( cell.name = colnames(dat@data), dat@meta.cell)
					namecol = 'cell.name'
				}
			}
			else {
				namecol = colnames(dat@meta.cell)[ok]
				namecol = namecol[1]
			}
			
			namerow = NULL
			if (nrow(dat@meta.gene)==0) {
				dat@meta.gene <- matrix(ncol=2, c(rownames(dat@data), rep( 0, nrow(dat@data)) ) )
				colnames(x@meta.gene) = c('Gene.Symbol', 'useless')
				rownames(x@meta.gene) = rownames(x@data)
				namerow = 'Gene.Symbol'
			}else {
				ok = which(lapply(colnames(dat@meta.gene) , function(x) { all.equal( as.character(as.vector(dat@meta.gene[,x])), rownames(dat@data)) == T } )==T)
				if ( length(ok) == 0) {
					if (all.equal( rownames(dat@meta.gene), rownames(dat@data)) == TRUE ){
						dat@meta.gene = cbind( gene.name = rownames(dat@data), dat@meta.gene)
						namerow = 'gene.name'
					}
				}
				else {
					namerow = colnames(dat@meta.gene)[ok]
					namerow = make.names(namerow[1])
				}
			}
			#storage.mode(dat@data) <- 'numeric'
			if ( nrow(dat@userGroups) == nrow(dat@meta.cell)){
				samples <- data.frame(cbind(dat@meta.cell, dat@userGroups))
			}else {
				samples <- data.frame(dat@meta.cell)
			}
			samples[,namecol] <- make.names(samples[,namecol])
			
			ret <- BioData$new( Matrix(dat@data), Samples=samples,
					annotation= dat@meta.gene, name= 'from.cellexalvr', namecol= namecol, namerow=namerow, outpath='./' )
			ret$usedObj <- dat@usedObj
			ret
		}
		
)

#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @param dat a seurat object
#' @title description of function as_BioData
#' @export 
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
						dat@hvg.info <- data.frame( 'Gene.Symbol' = rownames(dat@data),'useless' = rep( 0, nrow(dat@data) ) )
						rownames(dat@hvg.info) = rownames(dat@data)
						namerow = 'Gene.Symbol'
					}else {
						ok = which(lapply(colnames(dat@hvg.info) , function(x) { all.equal( as.character(as.vector(dat@hvg.info[,x])), rownames(dat@data)) == T } )==T)
						if ( length(ok) == 0) {
							## add the missing info
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
					#storage.mode(dat@data) <- 'numeric'
			
					dataa <- as.matrix(dat@data[1:2,1:2])
					samples <- data.frame(as.matrix(dat@meta.data))
					
					samples[,namecol] <- make.names(samples[,namecol])
					ret <- BioData$new( dataa, Samples=samples[1:2,], annotation =dat@hvg.info[1:2,], name= 'from.cellexalvr', namecol= namecol, namerow=namerow, outpath='./' )
					
					ret$dat <- Matrix(dat@data)
					ret$samples <- samples
					ret$annotation <- dat@hvg.info
					
					if ( ! is.null(dat@scale.data)){
						ret$zscored <- Matrix(as.matrix(dat@scale.data))
					}
					m <- match(colnames(dat@data), colnames(dat@raw.data))
					ret$raw <- Matrix(as.matrix(dat@raw.data))[,m]
					## now I need to manually remove unreliable data from the zscored information (set to -20)
#					if ( FALSE ) { # slow and not correct any more!
#						if ( ! is.null( ret$zscored )) {
#							d <- function( i, obj) {
#								drop <- which(obj$raw[,i] == 0)
#								if ( length(drop) > 0 ){
#									print ( paste( 'line', i, 'drop', length(drop)) )
#									ret$zscored[drop,i] = -20
#								}
#								0
#							}
#							lapply(  1:ncol(ret$dat), d, ret)
#						}
#					}
					ret$zscored= NULL
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


SQLite_2_matrix <- function ( fname, useS=NULL, useG=NULL, cells=  list( 'table' = 'samples', 'rev' = 'sample_id', 'name' = 'sname')) {
	dbh <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=fname )
	
	ncol=as.numeric( fetch_first( RSQLite::dbSendQuery(dbh,  paste("select max(id) from", cells$table ) ) ) )
	if ( ! is.null(useS) ) {
		ncol = length(useS)
	}
	nrow=as.numeric( fetch_first( RSQLite::dbSendQuery(dbh,  "select max(id) from genes" ) ) )
	if ( ! is.null(useG) ) {
		ncol = length(useG)
	}
	print ( paste("I create a", ncol,"columns and",nrow,"rows wide matrix"))
	
	q = paste( "select gene_id,", cells$rev,", value from datavalues")
	
	sth <- RSQLite::dbSendQuery(dbh, q )
	t <- RSQLite::dbFetch(sth)
	ret <- sparseMatrix( i=t[,1], j=t[,2], x=t[,3], dims=c(max(t[,1]),max(t[,2])))
	## gives the same result as the old function! EXTREMELY much faster!
	if(is.null(useS)){
		useS = 1:ncol(ret)
	}
	id = 1
	if (  is.null(useG) ){
		useG = 1:nrow(ret)
	}
	
	
	print ( "adding gene and sample names")
	RSQLite::dbClearResult(sth)
	rm(sth)
	q <- "select gname from genes order by id"
	sth <- RSQLite::dbSendQuery(dbh,q)
	rownames(ret) <- as.character(t(RSQLite::dbFetch(sth)))[1:nrow(ret)]
	RSQLite::dbClearResult(sth)
	rm(sth)
	
	q <- paste("select",cells$name,"from", cells$table, 'order by id' )

	sth <- RSQLite::dbSendQuery(dbh, q )
	colnames(ret) <- as.character(t(RSQLite::dbFetch(sth)))[1:ncol(ret)]
	RSQLite::dbClearResult(sth)
	RSQLite::dbDisconnect(dbh)
	print ( "done creating the matrix")
	ret
}

#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @param dat a sqlite database as created by the cellexalvrR::export2cellexalvr() function
#' @title description of function as_BioData
#' @export 
setMethod('as_BioData', signature = c ('character'),
		definition = function ( dat, minUMI=100, minGexpr=NULL ) {
			path = dat
			if ( ! file.exists( path ) ) {
				stop( "please - I need an existsing sqlite db at start up!")
			}
			print ( "Obtaining the sample UMI summary")
			
			sample_UMIs = tryCatch({ SQLite_SampleSummary( dat ) }, error = function(e) { 
						SQLite_SampleSummary( dat, list( 'table' = 'cells', name='cname', rev='cell_id') )
					} )
			useS =  sample_UMIs$sample_id[which(sample_UMIs$count > minUMI ) ]
			
			useG = NULL
			gene_UMI <- SQLite_ExpressionSummary( dat )
			if ( ! is.null(minGexpr)) {
				gene_UMI <- gene_UMI[which(gene_UMI$'count(value)' > minGexpr )]
			}
			useG = length(gene_UMI$gene_id)
			
			mat = tryCatch({ SQLite_2_matrix ( dat, useS = useS, useG = NULL ) } ,
					 error = function(e) {
						 SQLite_2_matrix ( dat, useS = useS, useG = NULL, list( 'table' = 'cells', name='cname', rev='cell_id') )
					 }
			 )
			
			print ( "creating the BioData object" )
			bad <- which( apply(mat,1, function(x) { length(which( x != 0)) } ) == 0)
			if ( length(bad) > 0 ) {
				mat <- mat[ -bad ,]
			}
			#mat$'ensembl_id' <- rownames(mat)
			gene_UMI <- gene_UMI[ match(rownames(mat),gene_UMI$gname ),]
			colnames(mat) <- make.names(colnames(mat))
			ret <- BioData$new( mat, annotation=gene_UMI ,
					Samples= data.frame( 'cell_id' = colnames(mat) ), namecol='cell_id', namerow= 'gname', name='from.cellranger' )
			class(ret) <- c( 'SingleCells', class(ret))
			changeNames(ret, 'row', 'gname')
			changeNames(ret, 'col', 'cell_id')
			ret
		}
)