fit_4_rf <- function ( x, copy= TRUE ) {
	if ( copy ) {
		x <- x$clone()
	}
	x$raw = NULL
	x$zscored = NULL
	m <- min(x$dat)
	if ( m == -1 ) {
		lapply( 1:ncol(x$dat), function(i) { 
					to0 <- which(x$dat[,i] == -1 ) 
					if ( length(to0) > 0 ){
						x$dat[to0,i] = 0
					}
					NULL
					 } 
			 )
	}
	invisible(x)
}

