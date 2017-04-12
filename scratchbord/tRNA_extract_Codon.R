stat_table <- function ( l ) {
	ret <- NULL
	for ( name in names(l)){
		if ( class(l[[name]]) == 'matrix'){
			ret <- rbind( ret, cbind( rep(name, nrow(l[[name]])), l[[name]]))
		}     
	}     
	ret   
}

