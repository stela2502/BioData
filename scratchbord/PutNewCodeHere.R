auto_reorder_grouping <- function( x, group ) {
	
	colapsed = collaps( x, groupCol= group , copy=T, by='sum')
	class(colapsed) = c('BioData', 'R6')
	colapsed$raw= colapsed$dat
	normalize(colapsed, force=T)
	
	clusters( colapsed, groups.n= 2 )
	new_order <- colapsed$usedObj$hc$labels[colapsed$usedObj$hc$order]
	browser()
	reorder_grouping( x, group= group, new_order= new_order )
	new_order
}
