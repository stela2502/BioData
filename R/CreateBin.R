#' @name CreateBin
#' @aliases CreateBin,BioData-method
#' @rdname CreateBin-methods
#' @docType methods
#' @description Bin the UMI data into 13 bins for plotting and define a blue <- red color gradient
#' @param x the BioData object
#' @param group the group (defaule = 'nUMI'
#' @param where in the samples (sample) or annotation (gene) data frame (sample)
#' @param colFun colory function default=  bluered
#' @title Create a binned annotation column from numeric data
#' @export 
setGeneric('CreateBin', ## Name
	function (x, group = 'nUMI', where='sample', colFun =  bluered ) { 
		standardGeneric('CreateBin')
	}
)

setMethod('CreateBin', signature = c ('BioData'),
	definition = function (x, group = 'nUMI', where='sample', colFun =  bluered ) {
	if ( where == 'sample' ){
		n <-as.numeric(as.vector(x$samples[,group] ))
	}else if ( where == 'gene' ) {
		n <-as.numeric(as.vector(x$annotation[,group] ))
	}else {
		stop(paste("Sorry where =",where,"is not supported (only sample and gene)") )
	}
	m <- min( n )
	brks= c( (m-.1),m ,as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
	brks = unique(as.numeric(sprintf("%2.6e", brks)))
	d  <- factor(brks [cut( n, breaks= brks)], levels=brks)
	newN = paste(sep="_",group, 'bins')
	if ( where == 'sample' ){
		x$samples[, newN] <- d
	}else {
		x$annotation[, newN] <- d
	}
	# defined the color
	x$usedObj$colorRange[[newN]] <- colFun( 13 )
	invisible(x)
} )


setMethod('CreateBin', signature = c ('BioData'),
		definition = function (x, group = 'nUMI', where='sample', colFun =  bluered ) {
			if ( where == 'sample' ){
				n <-as.numeric(as.vector(x$samples[,group] ))
			}else if ( where == 'gene' ) {
				n <-as.numeric(as.vector(x$annotation[,group] ))
			}else {
				stop(paste("Sorry where =",where,"is not supported (only sample and gene)") )
			}
			m <- min( n )
			brks= c( (m-.1),m ,as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
			brks = unique(as.numeric(sprintf("%2.6e", brks)))
			d  <- factor(brks [cut( n, breaks= brks)], levels=brks)
			newN = paste(sep="_",group, 'bins')
			if ( where == 'sample' ){
				x$samples[, newN] <- d
			}else {
				x$annotation[, newN] <- d
			}
			# defined the color
			x$usedObj$colorRange[[newN]] <- colFun( 13 )
			invisible(x)
		} )