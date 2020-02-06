readLoomFile <- function( x ) {
	require( 'hdf5r' )
	ret = as_BioData()
	
	x= loomR::connect(filename = x, mode = "r")
	
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
		dat <- x[['matrix']][, ]
	}
	# x will be an S3 matrix if X was scaled, otherwise will be a dgCMatrix
	if (is.matrix(x = dat)) {
		## the loom files seam to store all data as matrix and not as sparse matrix?!
		dat = Matrix::Matrix(dat, sparse=T)
	}
	print ( "reading feature data")
	ret$annotation = H5Anno2df(x, 'row_attrs', 'gene_names', onlyStrings=TRUE ) #function definition in file 'as_cellexalvrR.R'
	dat = Matrix::t(dat)
	rownames(dat) = rownames( meta.features)
	colnames(dat) = rownames(obs)
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
}
