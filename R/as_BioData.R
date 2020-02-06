#' @name as_BioData
#' @docType methods
#' @description create a BioData from an other object
#' @param ... source object specififc options
#' @title create a BioData from an other object
#' @export 
if ( ! isGeneric('as_BioData') ){ methods::setGeneric('as_BioData', ## Name
	function ( dat, ... ) { 
		standardGeneric('as_BioData') ## der Aufruf von standardGeneric sorgt für das_BioData Dispatching
	}
)
}else {
	print ("Onload warn generic function 'as_BioData' already defined - no overloading here!")
}

#' @describeIn as_BioData create a BioData::R6 object from a Rsubread result list
#' @docType methods
#' @param dat a list coming for Rsubread scanning of NGS data
#' @title create a BioData object from a Rsubread result list
#' @export as_BioData
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
				dat= Matrix::Matrix(dat$counts), 
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

setMethod("as_BioData","missing",function(x) {
			m = matrix(1,ncol=10, nrow=10)
			rownames(m) =paste( 'gene', 1:10)
			colnames(m) =paste( 'cell', 1:10)
			ret = as_BioData(m)
			ret$name="FakeObject"
			ret
		})

#' @describeIn as_BioData convert a matrix to a BioData::R6 object
#' @docType methods
#' @param dat a matrix object with cells in columns and genes in rows
#' @title convert a matrix into BioData
#' @export as_BioData
setMethod('as_BioData', signature = c ('matrix'),
		definition = function ( dat ) {
				dat <- Matrix::Matrix(dat)
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

#' @describeIn as_BioData convert a data.frame to a BioData::R6 object
#' @docType methods
#' @param dat a cellexlvrR object
#' @title convert a data.frame to BioData
#' @export as_BioData
setMethod('as_BioData', signature = c ('data.frame'),
		definition = function ( dat ) {
			snames <- colnames(dat)
			#dat$GeneID = rownames(dat)
			ret <- BioData$new( 
					dat=  Matrix::Matrix(as.matrix(dat)), 
					Samples = data.frame('sampleName' = snames ), 
					annotation= data.frame( 'GeneID' = rownames(dat)),
					namecol= 'sampleName', 
					namerow= 'GeneID',
					outpath= ''
			)
			ret
		}
)

#' @describeIn as_BioData convert a cellexalvrR object to a SingleCells::BioData::R6 object
#' @docType methods
#' @param dat a cellexlvrR object
#' @title convert a cellexalvrR object to SingleCells::BioData::R6
#' @export as_BioData
setMethod('as_BioData', signature = c ('cellexalvrR'),
		definition = function ( dat ) {
			#dat <- cellexalvr::renew(dat)
			#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
			m = as_BioData()

			ret$dat = dat@data
			ret$samples=cbind( 'cell.name' = colnames(ret$dat), dat@userGroups )
			ret$usedObj <- dat@usedObj
			ret$usedObj$Seurat.meta.cell = dat@meta.cell
			class(ret) = c( 'SingleCells', 'BioData', 'R6')

			return( ret)

			

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
			browser()
			if ( nrow(dat@userGroups) == nrow(dat@meta.cell)){
				samples <- data.frame(cbind(dat@meta.cell, dat@userGroups))
			}else {
				samples <- data.frame(dat@meta.cell)
			}
			samples[,namecol] <- make.names(samples[,namecol])
			
			ret <- BioData$new( Matrix::Matrix(dat@data), Samples=samples,
					annotation= dat@meta.gene, name= 'from.cellexalvr', namecol= namecol, namerow=namerow, outpath='./' )
			ret$usedObj <- dat@usedObj
			class(ret) = c( 'SingleCells', 'BioData', 'R6')
			ret
		}
		
)


#' @describeIn as_BioData convert a seurat object into a SingleCells::BioData::R6 object
#' @docType methods
#' @param dat a seurat object
#' @title convert a seurat object into BioData
#' @export as_BioData
setMethod('as_BioData', signature = c ('Seurat'),
				definition = function ( dat ) {
					
 					m = as_BioData() 
					ret$dat = dat@assays$RNA@data
					ret$raw = Seurat::GetAssayData(t, slot = "counts")
					nCell = FastWilcoxTest::ColNotZero( Matrix::t(ret$dat) )
					varG = rep( FALSE, nrow(ret$dat))
					varG[ match(  Seurat::VariableFeatures(dat), rownames(ret))] = TRUE
					ret$annotation= data.frame( GeneID = rownames(ret$dat), 'nCell' = nCell, varGene = varG )
					ret$samples = dat@meta.data
					ret$samples$nUMI = ret$samples$nFeature_RNA
					CreateBin( ret )
					ret$name = dat@project.name
					ret$logged = TRUE
					ret$snorm = TRUE
					class(ret) = c( 'SingleCells', 'BioData', 'R6')
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
		warning ("query was not finished")
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
	q = NULL
	if (  ! is.null(useS) ) {
		q = paste( "select gene_id,", cells$rev,", value from datavalues where ",cells$rev, " IN (", paste( useS, collapse = ", ") , ")")
	}else {
		q = paste( "select gene_id,", cells$rev,", value from datavalues")
	}
	
	sth <- RSQLite::dbSendQuery(dbh, q )
	t <- RSQLite::dbFetch(sth)
	ret <- Matrix::sparseMatrix( i=t[,1], j=t[,2], x=t[,3])
	
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
	#browser()
	if (  ! is.null(useS) ) {
		q <- paste("select",cells$name,"from", cells$table, " where id IN (", paste( useS, collapse= ", ") , ") ", 'order by id' )	 
	}else {
		q <- paste("select",cells$name,"from", cells$table,'order by id' ) 
	}
	sth <- RSQLite::dbSendQuery(dbh, q )
	t <- as.character(t(RSQLite::dbFetch(sth)))
	colnames(ret) <- t
	RSQLite::dbClearResult(sth)
	RSQLite::dbDisconnect(dbh)
	print ( "done creating the matrix")
	ret
}

load_database <- function( dat, minUMI=100, minGexpr=NULL ) {
	
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
				gene_UMI <- gene_UMI[which(gene_UMI$'count(value)' > minGexpr ),]
			}
			useG = length(gene_UMI$gene_id)
			mat = tryCatch({ SQLite_2_matrix ( dat, useS = useS, useG = NULL , list( 'table' = 'samples', name='sname', rev='sample_id')) } ,
					 error = function(e) {
						 SQLite_2_matrix ( dat, useS = useS, useG = NULL, list( 'table' = 'cells', name='cname', rev='cell_id') )
					 }
			 )
			
			print ( "creating the BioData object" )
			bad = which( Matrix::rowSums( mat ) == 0 )
			#bad <- which( apply(mat,1, function(x) { length(which( x != 0)) } ) == 0)
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

#' Supports sqlite3 databases as used in the cellexalVR application, 
#' 10X output files and kallisto out files.
#' 
#' @describeIn as_BioData character path or sqlite3 database name
#' @docType methods
#' @param dat a sqlite database as created by the cellexalvrR::export2cellexalvr() function
#' @param minUMI the minimum UMI count in order to read a cell (default=100)
#' @param minGexpr the minimum UMI count for a gene to be included (default NULL - translates into 1)
#' @title convert a cellRanger output directory or a Chromium_SingleCell_Perl output database into a BioData object
#' @export as_BioData
setMethod('as_BioData', signature = c ('character'),
		definition = function ( dat, minUMI=100, minGexpr=NULL ) {
			if ( dir.exists(dat) ){
				message("Process 10x output")
				return( load_10x(dat, minUMI=minUMI, minGexpr=minGexpr ) ) 
			}else { ## should be a sqlite databse then
				if ( length(grep( "loom$", dat))> 0 ){
					## readLoomFile
				}
				message("Load sqlite database")
				return( load_database(dat, minUMI=minUMI, minGexpr=minGexpr ) ) 
			}
			
		}
)

load_10x <- function ( dat, minUMI=100, minGexpr=NULL ) {
	err = NULL
	mtxFile = 'matrix.mtx'
	barFile = 'barcodes.tsv'
	geneFile = 'features.tsv'
	
	for( f in c('barcodes.tsv','matrix.mtx') ) {
		if ( ! file.exists(file.path(dat,f)) ) {
			if ( file.exists(file.path(dat,paste(f,sep='.', 'gz')))) {
				system( paste('gunzip ', file.path(dat,paste(f,sep='.', 'gz'))) )
			}else {
				err = c( err, paste("file", file.path(dat,f), "does not exists" ) )
			}
		}
	}
	OK = 0
	for( f in c('genes.tsv','features.tsv') ) {
		if ( file.exists(file.path(dat,f)) ) {
			OK = 1
			geneFile = f
		}
		if ( file.exists(file.path(dat,paste(f,sep='.', 'gz')))) {
			system( paste('gunzip ', file.path(dat,paste(f,sep='.', 'gz'))) )
			OK = 1
			geneFile = f
		}	
	}
	if ( OK == 0 ) {
		## check if the output is kallisto output
		OK = 1
		err2 = ""
		for ( f in c( 'genes.genes.txt', 'genes.barcodes.txt','genes.mtx')){
			if ( ! file.exists(file.path(dat,f)) ) {
				if ( file.exists(file.path(dat,paste(f,sep='.', 'gz')))) {
					system( paste('gunzip ', file.path(dat,paste(f,sep='.', 'gz'))) )
				}else {
					err2 = c( err2, paste("file", file.path(dat,f), "does not exists" ) )
					OK = 0
				}
			}
		}
		if ( OK ) {
			## load kallisto data
			mtxFile  = 'genes.mtx'
			barFile  = 'genes.barcodes.txt'
			geneFile = 'genes.genes.txt'
			err= NULL;
		}else {
			err = c( err, paste("file", file.path(dat,'features.tsv'), "does not exists" , err2) )
		}
	}
	
	if ( ! is.null(err)) {
		stop( paste(sep="\n", err))
	}
	head <- readLines( file.path(dat, mtxFile), 3)
	
	head =  unlist( stringr::str_split( head[3], ' ') )
	
	d = utils::read.table( file.path(dat,mtxFile) , comment.char="%", skip =1, header=F )
	
	ma = Matrix::sparseMatrix( i=d[,1], j=d[,2], x=d[,3] )
	message(paste("Matrix created with", ncol(ma) ,"cells and ", nrow(ma), "genes") )
	
	d = scan( file.path(dat,barFile) , what='character')
	if ( nrow(ma) == length(d) ) {
		ma = Matrix::t(ma)
	}
	colnames(ma) = d[1:ncol(ma)]
	
	rm(d)
	d= NULL
	d = utils::read.table( file.path(dat,geneFile), header=F )

	rownames(ma) = d[1:nrow(ma),1]
	
	message ( "Create BioData object" )
	ret = as_BioData(matrix( 1, ncol=2, nrow=2, dimnames= list( c('A','B'), c('C','D')) ))
	ret$dat = ma
	ret$samples = data.frame( 'samples' = colnames(ma), 'nUMI' = Matrix::colSums(ma) )
	if ( ncol(d) > 1 ){
		ret$annotation= data.frame( GeneID= d[1:nrow(ma),1], 'gene_name' = d[1:nrow(ma),2], 'nUMI' = Matrix::rowSums(ma))
		changeNames(ret, 'row', 'gene_name')
	}
	else {
		ret$annotation= data.frame( GeneID= d[1:nrow(ma),1], 'nUMI' = Matrix::rowSums(ma))
	}
	#ret$rownamescol = 'gene_name'
	
	reduceTo(ret, what='col', to=colnames(ret$dat)[which(ret$samples$nUMI >=  minUMI )], copy=F)
	if ( ! is.null(minGexpr) ) {
		reduceTo(ret, what='row', to=rownames(ret$dat)[which(ret$annotation$nUMI >=  minGexpr )], copy=F)
	}
	
	rm(ma)
	rm(d)
	gc()
	message( paste("Final object with",ncol(ret$dat),"cells and", nrow(ret$dat),"genes") )
	return( ret )
	
}
