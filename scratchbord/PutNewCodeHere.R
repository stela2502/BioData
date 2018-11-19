cells_overlap_group <- function( x, A, colName, from='data') {
	ret = data.frame(lapply( levels(GOI1.2$samples[,colName] ),
		function(i){
			cells_overlap_cells( GOI1.2, A, which( GOI1.2$samples[,colName] == i), from="Expression PCA") } 
		))
	colnames(ret) = levels(GOI1.2$samples[,colName] )
	rownames(ret) = c('max', 'min', 'mean', 'median')
	if ( all( is.na( as.vector(t(ret['mean',])) ) ) ){
		return ( "noMatch" )
	}
	m = min(as.vector(t(ret['mean',])), na.rm=T)
	colnames(ret)[which(ret['mean',] == m)]
}


cells_overlap_cells <- function( x, A, B, from='data' ) {
	if ( length(A) < 3 | length(B) < 3) {
		return (c(NA, NA, NA, NA))
	}
	locA = cell_locations(x, A, from=from)
	locB = cell_locations(x, B, from=from)
	dist = abs(locA[,'mean'] - locB[,'mean'] )
	ok = c( max( dist ), min(dist), mean(dist), median(dist))
	ok
}

cell_locations <- function(x, cells, from='data' ) {
	## lets assume we use the data only
	if ( from == 'data' ) {
		ret = t( apply( x$data()[,cells],1,function(x) { c(mean(x), sd(x)) }))
		colnames(ret) = c('mean', 'sd')
	}else if ( match( from, names(x$usedObj$'MDS_PCA100'))){
		ret = t( apply( x$usedObj$'MDS_PCA100'[[from]][cells,],2,function(x) { c(mean(x), sd(x)) }))
		colnames(ret) = c('mean', 'sd')
	}
	else {
		stop("Sorry only data is supported at the time" )
	}
	ret
}