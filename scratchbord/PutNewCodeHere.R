
force.numeric = function ( x ) {
	for ( i in 1: ncol(x$data) ) {
		if ( ! is.numeric(x$data[,i]) ) {
			x$data[,i] <- as.numeric(x$data[,i])
		}
	}
}