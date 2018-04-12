add_to_stat <- function( x, stat, name ) {
	if ( length( match('p.adj BF',colnames(stat) )) == 0 & length( match('hurdle',colnames(stat) ))) {
		stat = cbind( stat, 'p.adj BF' = p.adjust(stat[,'hurdle'], method='bonferroni') )
	}
	if ( ! is.na( match( name, names(x$stats)))){
		x$stats[[ match( name, names(x$stats)) ]] <- stat
	}else {
		x$stats[[ length( x$stats ) +1 ]] <- stat
		names(x$stats)[length(x$stats) ] <- name
	}
	invisible(x)
}
