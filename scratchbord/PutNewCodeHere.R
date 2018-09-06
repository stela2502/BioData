group_intersect_order <- function( x, group, otherGroup ) {
	r = lapply( levels(x$samples[, group]), 
			function(n) { 
				t = table(x$samples[, otherGroup][which(x$samples[, group] == n )])
				c(names(t[which( t == max(t))]),n) 
			} 
	)
	r= t(data.frame(r))
	rownames(r) = NULL
	new_order = order(factor( r[,1], levels( x$samples[, otherGroup] )))
	reorder_grouping(x, group= group,  new_order=as.vector(r[new_order,2]),  what='row')
	invisible(x)
} 
