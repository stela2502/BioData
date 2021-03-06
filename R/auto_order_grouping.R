#' @name auto_order_grouping
#' @aliases auto_order_grouping,BioData-method
#' @rdname auto_order_grouping-methods
#' @docType methods
#' @description This function uses the rfCluster_col function to order median collapsed data values.
#' This RF order is then re-ordered twice using hclust 'ward.D2' and the final ordering is forced onto the original grouping.
#' @param x The BioData object
#' @param group the sample grouping that should become reordered (min 15 groups)
#' @param settings a list of slurm settings default list( 'A' = 'lsens2017-3-2' , 't'= '01:00:00', 'p' = 'dell', 'n' = 1, 'N' = 1)
#' @param k how many groups to put into one super group (default 5)
#' @param colapseBy which method to collapse the data by (default 'median') see function collaps 
#' @title Auto re-order a (RF) grouping
#' @export 
setGeneric('auto_order_grouping', ## Name
	function ( x, group,  settings=list( 'A' = 'lsens2017-3-2', 't' = '01:00:00', p='dell', 'n'=1, 'N'=1) , k=5,
			colapseBy = 'median') { 
		standardGeneric('auto_order_grouping')
	}
)

setMethod('auto_order_grouping', signature = c ('BioData'),
	definition = function ( x, group,  settings=list( 'A' = 'lsens2018-3-3', 't' = '01:00:00', p='dell', 'n'=1, 'N'=1), k=5, 
			colapseBy = 'median' ) {
	
	k = ceiling(length( unique(x$samples[,group])) / k )
	if ( k < 3) {
		stop( "too view groups to try an automatic grouping here" )
	}
	print ("Collapsing the data")
	x_collapsed = collaps( x, by= colapseBy , groupCol= group, copy= TRUE )
	x_collapsed$name='AutoOrder1'
	
	tab <- x_collapsed$data()
	p <- apply(tab, 2, var)
	problem = which(p == 0 | is.na(p) )
	if ( length(problem) > 0 ){
		stop(paste("The group(s)", paste( collapse=', ', names(problem)), "have a var of 0 or NA in the collapsed dataset - please add more genes and re-run or cange the colapseBy option to 'sum'") )
	}
	
	print ( "Starting RF clustering process")
	rfCluster_col(x_collapsed, k=k, slice=20, subset = ncol(x_collapsed$dat), name='autoorder1', settings=settings )
	exp_group = paste('RFgrouping autoorder1 .*', ' n=', k, sep="")

	while ( length( grep ( exp_group, colnames( x_collapsed$samples )) ) == 0){
		print("wait for the RF process to finish" )
		Sys.sleep( 30 )
		try( {rfCluster_col(x_collapsed, k=k, slice=20, 
							subset = ncol(Stat_collapsed$dat), 
							name='autoorder1', 
							settings=settings )} , silent=T)	
	}
	exp_group = colnames( x_collapsed$samples )[grep ( exp_group, colnames( x_collapsed$samples ))[1]]
	x_collapsed2 <- collaps( x_collapsed, by='sum', groupCol= exp_group, copy= TRUE )
	print ( "creating hclust orders")
	tab <- x_collapsed2$data()
#	browser()
	colnames(tab) <- x_collapsed2$samples[,exp_group]
#	p <- apply(tab, 2, var)
#	problem = which(p == 0 | is.na(p) )
#	if ( length( problem ) > 0 ) {
#		tab = tab[,-problem]
#	}
	hc <- stats::hclust(as.dist( 1- stats::cor(tab, method='pearson') ),method = 'ward.D2' )
	order_l2 <- hc$order ## should be the most likely usable order here
#	if ( length( problem ) > 0 ) {
#		for ( i in problem ) { 
#			insert(order_l2, i, length(order_l2)+1) # push the empty groups to the end
#		}
#	}
	tab <- x_collapsed$data()
	colnames(tab) <- x_collapsed2$samples[,group]
	hc <- stats::hclust(as.dist( 1- stats::cor(tab, method='pearson') ),method = 'ward.D2' )
	order_l1 <- hc$order ## should be the most likely usable order here
	
	## now create the most likely OK new order :-D
	
	info <- split( as.vector(x_collapsed$samples[,group]) , x_collapsed$samples[,exp_group])
	
	neq_order <- NULL;
	for ( pos2 in order_l2 ) {
		neq_order <- c( neq_order, info[[pos2]][order(match(info[[pos2]], order_l1))])
	}
	print ( "apply new order and return" )
	reorder_grouping( x, group, neq_order )
	
	invisible(x)
} )
