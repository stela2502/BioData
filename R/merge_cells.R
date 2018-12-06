#' @name merge_cells
#' @aliases merge_cells,BioData-method
#' @rdname merge_cells-methods
#' @docType methods
#' @description merge cells to get summary values of median values or anything you want 
#' @param dataObj the BioData object
#' @param group the original grouping you want to merge on
#' @param mergeNcells merge mergeNcells into one summary cell (default =10)
#' @param gID merge only cells in a spicific ID list of this group; Say the group as entries 'A', 'B', and 'C', but you only want to merge cells on groupID 'A' and 'B' you need to specifiy ('A', 'B') here (default =NULL)
#' @param by  by collapsing method c('median','mean','sd','sum', or own function )
#' @param copy create a copy of the R6 object instead of changing the real thing (default=TRUE)
#' @title description of function merge_cells
#' @export 
setGeneric('merge_cells', ## Name
	function ( dataObj, group, mergeNcells=10, gID=NULL, by=c('median','mean','sd','sum', 'var' ), copy=TRUE ) { 
		standardGeneric('merge_cells')
	}
)

setMethod('merge_cells', signature = c ('BioData'),
	definition = function ( dataObj, group, mergeNcells=10, gID=NULL, by=c('median','mean','sd','sum', 'var' ), copy=TRUE ) {
	if ( copy ) {
		dataObj = dataObj$clone()
	}
	if ( is.function(by)){
		f <- by
		by = f  
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
		m = sample(m) ## a little bit more randomness is likely good
		if ( is.na(match( gname, gID) ) ){
			new_group[m] = i:(i+length(m) -1)
			i = i+length(m)
		}else {
			## OK here we need to merge the cells into summary samples by mergeNcells
			o <- split(m, ceiling(seq_along(m)/mergeNcells))
			for ( id in 1:length(o) ){
				new_group[o[[id]]] = i
				i = i+1
			}
		}
		## And now I can add this new column
		dataObj$samples$'MergeByThat' <- new_group
	}
	
	collaps(dataObj, groupCol='MergeByThat', by = by ) ## does not need to be copied as the copying has been done here.
	
} )
