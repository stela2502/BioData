#' @name readLoomFile
#' @aliases readLoomFile,BioData-method
#' @rdname readLoomFile-methods
#' @docType methods
#' @description read a loom file
#' @param x the loom file
#' @title description of function readLoomFile
#' @export 
if ( ! isGeneric('readLoomFile') ){setGeneric('readLoomFile', ## Name
	function ( x ) { 
		standardGeneric('readLoomFile')
	}
) }

setMethod('readLoomFile', signature = c ('character'),
	definition = function ( x ) {
	require( 'hdf5r' )
	ret = as_BioData()
	
	x= hdf5r::h5file(x)
	message(class(x))
	#x= loomR::connect(filename = x, mode = "r")
	
	as.sparse <- function(x, ...) {
		for (i in c('data', 'indices', 'indptr')) {
			if (!x$exists(name = i) || !is(object = x[[i]], class2 = 'H5D')) {
				stop("Invalid H5Group specification for a sparse matrix, missing dataset ", i)
			}
		}
		if ('h5sparse_shape' %in% hdf5r::h5attr_names(x = x)) {
			return(Matrix::sparseMatrix(
							i = x[['indices']][] + 1,
							p = x[['indptr']][],
							x = x[['data']][],
							dims = rev(x = hdf5r::h5attr(x = x, which = 'h5sparse_shape'))
					))
		}

		return(Matrix::sparseMatrix(
						i = x[['indices']][] + 1,
						p = x[['indptr']][],
						x = x[['data']][],
						sparse=T
				))
	}
	print ( "reading cell information" )
	ret$samples = obs = H5Anno2df(x,'col_attrs', 'cell_names', onlyStrings=TRUE ) #function definition in file 'as_cellexalvrR.R'
	
	print ( "reading data" )
	if (is(object = x[['matrix']], class2 = 'H5Group')) {
		dat <- as.sparse(x = x[['matrix']])
	} else {
		message( "Converting full matrix to sparse (SLOW)")
		obj = Dense2SparseHDF5::Dense2SparseHDF5$new('test')
		obj$file= x
		obj$toSparseVector()
		dat = obj$Matrix
		rm(obj)
		gc()
		message('finished')
	}
	# x will be an S3 matrix if X was scaled, otherwise will be a dgCMatrix

	print ( "reading feature data")
	ret$annotation = H5Anno2df(x, 'row_attrs', 'gene_names', onlyStrings=TRUE ) #function definition in file 'as_cellexalvrR.R'
	samples =  H5Anno2df(x, 'col_attrs',  onlyStrings=TRUE ) #function definition in file 'as_cellexalvrR.R'
	dat = Matrix::t(dat)
	browser()
	rownames(dat) = rownames( res$annotation )
	colnames(dat) = rownames(res$samples)
	ret$dat = dat
	rm( dat)

	## The drc's are hidden in the obj
	interest <- list( 
			'unknown' = c('_X', '_Y', '_Z'), 
			'tSNE' = c('_tSNE1', '_tSNE2', '_tSNE3'), 
			'PCA' = c('_PC1', '_PC2', '_PC3') 
	)
	print ("reading drc data")
	dr = lapply( interest, function(a) { 
				d=NULL
				if ( var(ret$samples[,a[1]]) + var (retsamples[,a[2]]) != 0 ){
					## the third column is not defined in the loom structures. Hence I simply do not check for it here
					d= cbind(as.numeric(as.vector(ret@usedObj$original_meta.cell[,a[1]])), as.numeric(as.vector(ret@usedObj$original_meta.cell[,a[2]])), rep(0,nrow(ret@meta.cell)) )
					colnames(d) = a
					rownames(d) = rownames(ret@meta.cell)
				}
				d
			} )
	
	ret$usedObj$MDS_PCA100 = dr
	
	ret$outpath = getwd()
	#print ( "Please take care colnames and rownames have not been set!" )
	invisible(ret)
} )


#' @name H5Anno2df
#' @aliases H5Anno2df,cellexalvrR-method
#' @rdname H5Anno2df-methods
#' @docType methods
#' @description  convert a H5 annotation (any name) table to a data table
#' @param x the H5 object
#' @param slotName the H5 entity tro convert to a data.frame
#' @param namecol the (optional) rownames column for the data
#' @param onlyStrings return only columns that not only contain numbers (default FALSE)
#' @title description of function H5Anno2df
#' @export 
setGeneric('H5Anno2df', ## Name
	function (x, slotName, namecol=NULL, onlyStrings=FALSE ) { ## Argumente der generischen Funktion
		standardGeneric('H5Anno2df') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('H5Anno2df', signature = c ('H5File'),
	definition = function (x, slotName, namecol=NULL, onlyStrings=FALSE ) {
  		obs = data.frame(lapply(names(x[[slotName]]), function(n) { x[[paste(sep="/",slotName,n)]][] } ))
  		colnames( obs ) = names(x[[slotName]])
  		col_uniq= NULL
  		for( n in colnames(obs) ) {
	  		if ( all(obs[,n] =="") ){
  				obs[,n] = NULL
  			}else {
  				col_uniq = c( col_uniq, length(unique(obs[,n]))) 
  			}
  		}
  		names(col_uniq) = colnames( obs )
  		## most likely cell names column
  		if ( ! is.na(match(namecol, colnames(obs)) )) {
			rownames(obs) =  make.unique ( #function definition in file 'as_cellexalvrR.R'
				as.vector(obs[, namecol]) )
  		}else {
  			## now I need to check for strings...
  			OK = unlist(lapply( colnames(obs) , function(id) {
  				a= which( is.na(as.numeric(as.vector(obs[,id])))==T) ## Strings only
  				if ( length(a) > 0) {
  					length(unique(as.vector(obs[a, id])))
  				}else {
  						0
  				}
  			}))
  			names(OK) = colnames(obs)
  			# if ( slotName == 'row_attrs'){
  			# 	browser()
  			# }
  			rownames(obs) =  make.unique ( #function definition in file 'as_cellexalvrR.R'
  				as.vector(obs[, names(OK)[which(OK == max(OK))[1]]]) )
  		}
  		if ( onlyStrings ) {
  			for( i in 1:length(col_uniq) ) {
  				if ( col_uniq[i] == 0 ){ # all entries convertable to numeric
  					obs[,i] = NULL
  				}
  			}
  		}

  		obs
  	} )