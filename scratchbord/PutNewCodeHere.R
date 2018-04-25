merge_cells <- function( dataObj, group, mergeNcells=10, gID=NULL, by=c('median','mean','sd','sum', 'var' ), copy=TRUE ) {
	if ( copy ) {
		dataObj = dataObj$clone()
	}
	if ( is.function(by)){
		f <- by
		by = "user_function"
	}else {
		switch( by,
				median = f<- function (x ) { median(x) },
				mean = f <- function(x) { mean(x) },
				sd = f <- function(x) { sd(x) },
				sum = f <-function(x) { sum(x)},
				var = f <- function(x) { var(x) }
		);
	}
	if ( is.null(f) ) {
		stop("Please set what to one of 'median','mean','sd','sum'" )
	}
	
	new_samples <- NULL
	
	if ( is.null(gID) ) {
		gID = unique(as.vector(dataObj$samples[,group]))
	}
	
	## how many new samples do I need?
	gIDs = unique(as.vector(dataObj$samples[,group]))
	new_group = rep('problem', nrow(dataObj$samples))
	i =1
	for( gname in gIDs ) {
		m <- match( dataObj$samples[,group], gname )
		m <- which(is.na(m) == F)
		
		if ( is.na(match( gname, gID) ) ){
			new_group[m] = 1:(i+length(m))
			i = i+length(m)
		}else {
			## OK here we need to merge the cells into summary samples by mergeNcells
			o <- split(m, ceiling(seq_along(m)/mergeNcells))
			for ( i in 1:length(o) ){
				new_group[m[[o]]] = i
				i = i+1
			}
		}
		## And now I can add this new column
		dataObj$samples$'MergeByThat' <- new_group
	}
	
	collaps(dataObj, groupCol='MergeByThat', by = by ) ## does not need to be copied as the copying has been done here.
	
}