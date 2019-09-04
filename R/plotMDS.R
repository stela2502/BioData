#' @name plotMDS
#' @aliases plotMDS,BioData-method
#' @rdname plotMDS-methods
#' @docType methods
#' @description Create 2D plts of any MDS and a list of feature plots
#' @param obj The BioData object
#' @param mds The name of the MDS to plot
#' @param group Which sample grouping(s) to plot onm the MDS  
#' @param genes Which genes to color the MDS with
#' @param x MDS dimension on the x axis default=1
#' @param y MDS dimension on the y axis default=2
#' @param green show the gene expression lost during normalization in green default=F
#' @title description of function plotMDS
#' @export 
if ( ! isGeneric('plotMDS') ){setGeneric('plotMDS', ## Name
	function (obj, mds, group, genes, x=1, y=2, green=F ) { 
		standardGeneric('plotMDS')
	}
) }

setMethod('plotMDS', signature = c ('BioData'),
	definition = function (obj, mds, group, genes, x=1, y=2, green=F ) {
	mds.dat = NULL
	if ( ! is.null(obj$usedObj$MDS_PCA100_dim_2[[mds]] )){
		print( "Propper 2D dimension reduction found" )
		mds.dat = obj$usedObj$MDS_PCA100_dim_2[[mds]]
	}else if ( ! is.null(obj$usedObj$MDS_dim_2[[mds]] )) {
		print( "Propper 2D dimension reduction found" )
		mds.dat = obj$usedObj$MDS_dim_2[[mds]]
	}
	else if ( ! is.null(obj$usedObj$MDS[[mds]] )) {
		print( "3D dimension reduction object using dim 1 and 2" )
		mds.dat = obj$usedObj$MDS[[mds]]
	}
	else if ( ! is.null(obj$usedObj$MDS_PCA100[[mds]] )){
		print( "3D dimension reduction object using dim 1 and 2" )
		mds.dat = obj$usedObj$MDS_PCA100[[mds]]
	}
	fname <- function( parts ) {
		n = paste(parts, collapse="_")
		paste(stringr::str_replace_all( n, "\\s\\s+",'_'),'svg', sep='.')
	}
	
	if ( is.null(mds.dat)) {stop("MDS data not part of this object!") }
	for ( g in group ) {
		colors_4(obj, g) ## if the color is not already defined do it here
		svg( file.path( obj$outpath, fname( c(obj$name, mds, x, y, g) )), width=10 , height=10)
		o = order(obj$samples[,g])
		plot ( mds.dat[o,x], mds.dat[o,y], col=obj$usedObj$colorRange[[g]][obj$samples[o,g]], xlab=paste('dim',x), ylab=paste('dim', y), pch=16 )
		dev.off()
	}
	m = min(obj$data())
	getF <- function(v) {
		l = levels(factor( v ))
		if ( length(which( l == 0)) == 0 ) {
			l= c(0,l)
		}
		if ( length(which( l == m)) == 0 ) {
			l= c(m,l)
		}
		factor( v, levels=l)
	}
	
	for ( g in genes ) {
		n <- as.numeric(obj$data()[g,] )
		col = NULL
		COLS = NULL
		d = NULL
		if ( m == -1 | m == -21){
			if ( length(table(n)) < 20 ) {
				## no need to break that into bins
				d = getF(n)
			}else {
				brks= c( (m-.1),m+1-.1 , as.vector(quantile(n[which(n > 0)],seq(0,1,by=0.1)) ) )
				brks = unique(as.numeric(sprintf("%2.6e", brks)))
				brks[3] = brks[3] - 1e-5
				brks[length(brks)] = brks[length(brks)]  + 0.1
				d  <- factor(brks [cut( n, breaks= brks)], levels=brks)
			}
			if ( green) {
				COLS = c('#006D2C', 'black', gplots::bluered(length(levels(d)) -2  ))
			}else{
				COLS = c('black', 'black', gplots::bluered(length(levels(d)) -2  ))
			}
			
			if ( length(COLS) == 3 ) {
				COLS[3] = 'blue'
			}
			if ( length(COLS) == 4 ) {
				COLS[4] = 'red'
			}
			col = COLS[d]
		}else {
			brks= c( (m-.1),m,as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
			brks = unique(as.numeric(sprintf("%2.6e", brks)))
			d  <- factor(brks [cut( n, breaks= brks)], levels=brks)
			COLS = c('black', gplots::bluered(length(brks) -1  ))
			if ( length(COLS) == 2 ) {## the second col is white - not good!
				COLS[2] = 'blue'
			}
			col = COLS[d]
		}
		svg( file.path( obj$outpath, fname( c(obj$name, mds, x, y, g) )), width=10, height=10)
		d = order(d)
		bad = grep('#FFFFFF', col)
		if ( length(bad) > 0 ){ col[bad] = '#FFDBDB' }
		plot ( mds.dat[d,x], mds.dat[d,y], col=col[d], xlab=paste('dim',x), ylab=paste('dim', y), pch=16, main=g )
		dev.off()
	}
} )
