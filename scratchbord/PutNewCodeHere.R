
chromosomal_distance_from <- function(x, chr, pos, chrCol, posCol ){
	d = data.frame( id= 1:nrow(x$dat), chr=as.vector( x$annotation[,chrCol]), dist = abs( x$annotation[,posCol] - pos  ) )
	d=d[which(d[,2] == chr), ]
	d= d[order[d[,3]],]
	d
}