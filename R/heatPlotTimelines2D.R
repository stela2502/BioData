#' Plot summed up gene expression on top of a 2D dimension reduction a 
#' sample grouping and a timeline.
#' 
#' @name heatPlotTimelines2D
#' @aliases heatPlotTimelines2D,BioData-method
#' @rdname heatPlotTimelines2D-methods
#' @docType methods
#' @description Plot a summed up gene expression on top of a 2D dimension reduction and 
#' @param obj  TEXT MISSING
#' @param times  TEXT MISSING
#' @param mds  TEXT MISSING
#' @param grouping  TEXT MISSING
#' @param genes  TEXT MISSING
#' @param outpath  TEXT MISSING default=NULL
#' @param sumN  TEXT MISSING default=10
#' @param x  TEXT MISSING default=1
#' @param y  TEXT MISSING default=2
#' @title description of function heatPlotTimelines2D
#' @export 
if ( ! isGeneric('heatPlotTimelines2D') ){setGeneric('heatPlotTimelines2D', ## Name
	function ( obj, times, mds, grouping, genes, outpath=NULL, sumN=10, x=1, y=2  ) { 
		standardGeneric('heatPlotTimelines2D')
	}
) }

setMethod('heatPlotTimelines2D', signature = c ('BioData'),
	definition = function ( obj, times, mds, grouping, genes, outpath=NULL, sumN=10, x=1, y=2  ) {
	if ( ! is.null(obj$usedObj$MDS_PCA100_dim_2[[mds]] )){
		print( "Propper 2D dimension reduction found" )
		mds.dat = obj$usedObj$MDS_PCA100_dim_2[[mds]]
	}else if ( ! is.null(obj$usedObj$MDS_dim_2[[mds]] )) {
		print( "Propper 2D dimension reduction found" )
		mds.dat = obj$usedObj$MDS_dim_2[[mds]]
	}
	else if ( ! is.null(obj$usedObj$MDS[[mds]] )) {
		print( paste("3D dimension reduction object using dim",x, "and",y ) )
		mds.dat = obj$usedObj$MDS[[mds]]
	}
	else if ( ! is.null(obj$usedObj$MDS_PCA100[[mds]] )){
		print( paste("3D dimension reduction object using dim",x, "and",y ) )
		mds.dat = obj$usedObj$MDS_PCA100[[mds]]
	}
	if ( is.null(outpath) ) {
		outpath = obj$outpath
	}
	## probably not roll, but hop?!
	roll = function( x, smooth, type ) {
		func = NULL
		if ( type=='mean') {
			func = function(start, x, smooth ) { ret= mean( as.vector(x[(start-smooth):start])); if ( is.na(ret)){browser();}; ret  } 
		}else if ( type == 'table' ) {
			func = function( start, x, smooth) { length(names(table(as.vector(x[(start-smooth):start])))) }
		}else {
			stop ( "the roll function needs either 'mean' or 'table' as type!" )
		}
		#browser()
		c( rep(NA, smooth), unlist(lapply ( (smooth+1):length(x), func, x, smooth ) ) )
	}
	hop = function( x, smooth, type ) {
		func = NULL
		if ( type=='mean') {
			func = function(start, x, smooth ) { ret= mean( as.vector(x[start:(start+smooth-1)]), na.rm=T); rep(ret, smooth)  } 
		}else if ( type == 'table' ) {
			func = function( start, x, smooth) { ret = length(names(table(as.vector(x[start:(start+smooth-1)])))); rep( ret, smooth) }
		}else {
			stop ( "the roll function needs either 'mean' or 'table' as type!" )
		}
		#browser()
		ret= unlist(lapply ( seq( 1, length(x), smooth ), func, x, smooth ) )
		ret[1:length(x)]
	}
	
	dat = obj$dat
	dat@x[which(dat@x == -1)] = 0
	Matrix::drop0(dat)
	
	getData <- function( name )  {
		#browser()
		good = which(is.na( obj$samples[, paste(sep="_", name, 'Pseudotime')]) == F)# cells part of this timeline
		o = order(obj$samples[good, paste(sep="_", name, 'Pseudotime')])
		X = hop( dat[gname, good[o]], sumN, 'mean' )
		return( list( 
			X =X, 
			x= obj$samples[good[o], paste(sep="_", name, 'PseudotimeX')],
			y= obj$samples[good[o], paste(sep="_", name, 'PseudotimeY')]
		))

	}
	plotDat <- function(x, brks) {
		d= cut( x$X, brks)
		o = order(d)
		o = length(d):1
		#browser()
		points( as.vector(x$x[o]), as.vector(x$y[o]), col= 'black', pch=16, cex=2 )
		points( as.vector(x$x[o]), as.vector(x$y[o]), col= x$col[o], pch=16, cex=1.8 )
	}
	updateValues <- function( x, new, col) {
		start = 1
		for ( i in 1:length(x) ) {
			x[[i]]$X = new[start:(length(x[[i]]$X) + start -1 )]
			x[[i]]$col = col[start:(length(x[[i]]$X) + start -1 )]
			start = length(x[[i]]$X) + start
		}
		x
	}
	
	for ( gname in genes ) {
		ofile=paste( sep="_", obj$name,gname,mds)
		ofile= stringr::str_replace_all( ofile, "\\s+", '_' )
		ofile=paste( sep=".", ofile, 'pdf')
		message( paste( "create", ofile) )
		
		## get the data
		b=lapply( times, getData )
		names(b) = times
		all_values = unlist(lapply( b, function(x){x$X }))
		## z.score
		on = which( all_values != 0) # they should get colors!
		off = which( all_values == 0)
		all_values[on] = (all_values[on] - mean(all_values[on]))/ sd(all_values[on])
		all_values[which(all_values < -3)] = -3
		all_values[which(all_values > 3)] = 3
		all_values[off] = -4
		brks <- c(-30, unique(as.vector(c(0, stats::quantile(all_values,seq(0,1,by=1/100)),max(all_values)))))
		d = cut( all_values, brks)
		col = c('black', gplots::bluered(length(levels(d))-1) ) [d]
		
		b = updateValues(b , all_values, col)
		#brks = seq( min(all_values), max(all_values),  (max(all_values)- min(all_values))/50)
		#browser()
		#brks <- c(-30, unique(as.vector(c(0, stats::quantile(all_values,seq(0,1,by=1/100)),max(all_values)))))
		#browser()
		#brks[1] = brks[1] - 1e-4
		## plot the original grouping
		grDevices::pdf( file=file.path(outpath, ofile), width=20, height=20 )
		
		plot ( mds.dat[o,x], mds.dat[o,y], col=obj$usedObj$colorRange[[g]][obj$samples[o,g]], xlab=paste('dim',x), ylab=paste('dim', y) , pch=16, cex=3)
		lapply( b , plotDat, brks )
		dev.off()
	}
	invisible(obj)
} )
