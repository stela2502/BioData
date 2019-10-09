#' A simple pseudotime analysis based on a linear model, that can create the pseudotime from the two varibales 'a' and 'b'
#' 
#' fit = loess( b~a )
#' 
#' The fitted values will be used as pseudotime and the pseudotiume is plotted 
#' versus the two defininge variables to the file outpath/Pseudotime.png using the grouping color.
#' Pseudotime is depicted as bluered line indicating the orientation blue -> red.
#' 
#' The pseudotime will be added into the samples table and the correlation to 
#' the pseudotime will be added to the annotation table (both replacing old data!).
#' 
#' The function returns the top/bottom (n) correlating genes.
#' 
#' The rolling mean values that are summing up only expression values from one group 
#' are colored in the group colour all others are black. 
#'  
#' @name pseudotimeTest
#' @aliases pseudotimeTest,BioData-method
#' @rdname pseudotimeTest-methods
#' @docType methods
#' @description Calculate and plot pseudotime relevant information.
#' @param x a BioData object (preferably a SingleCells object)
#' @param a the first time defining vector (e.g. a MDS dimension)
#' @param b the second time defining vector (e.g. a MDS dimension)
#' @param grouping a samples grouping to color the Pseudotime plot.
#' @param outpath a outpath for the files  default=x$outpath
#' @param n how many genes to plot on both sides default=100
#' @param plotGenes additional default genes to plot default=NULL
#' @param smooth roling smoothing window size (default 100)
#' @param invert invert the time before usage ( default=FALSE)
#' @param cleanFolder remove all files from the outpout folder before creating new ones (default FALSE)
#' @param plotType ( 'png', 'png_high_res', 'pdf' )
#' @param summaryPlot the name of the summary plot file ( default NULL no summary plot)
#' @title description of function pseudotimeTest
#' @export 
if ( ! isGeneric('pseudotimeTest') ){setGeneric('pseudotimeTest', ## Name
			function ( x, a, b, grouping, outpath=NULL,  n=100, plotGenes=NULL, smooth=100, 
					invert=FALSE, cleanFolder=FALSE, plotType='png', summaryPlot=NULL ) { 
				standardGeneric('pseudotimeTest')
			}
	) }

setMethod('pseudotimeTest', signature = c ('BioData'),
		definition = function ( x, a, b, grouping, outpath=NULL,  n=100, plotGenes=NULL, 
				smooth = 100, invert=FALSE, cleanFolder=FALSE, plotType='png' , summaryPlot=NULL ) {
			
			openPlot <- function(fname) {
				if ( plotType == 'pdf' ) {
					grDevices::pdf( file=paste(fname,'pdf', sep="."), width=10, height=10 )
				}else if (plotType == 'png_high_res' ){
					grDevices::png ( file=paste(fname,'highRes','png', sep="."), width=1600, height=1600)
				}else {
					grDevices::png ( file=paste(fname,'png', sep="."), width=800, height=800)
				}
			}
			Ra = range(a)
			Rb = range(b)
			rev.order = FALSE
			if ( abs( Rb[1] - Rb[2]) > abs(Ra[1] - Ra[2]) ) {
				rev.order=TRUE
			}
			if ( is.null( outpath )) {
				outpath = x$outpath
			}
			if ( ! file.exists( outpath )) {
				dir.create( outpath )
			}else if ( cleanFolder) {
				unlink(file.path( outpath ,"*.png") )
				unlink(file.path( outpath ,"*.pdf") )
			}
			if ( is.null( names(a) )) {
				stop( "This function needs names on the a vector" )
			}
			#browser()
			pred2D <- function( a, b ) {
				fit = list()
				fit$a = a
				ls = loess( b ~ a )
				fit$b = predict( ls)
				fit$fitted.values = fit$b
				o = order(a)
				if ( invert ) {	o = rev(o)	}
				fit$time = FastWilcoxTest::euclidian_distances( a[o], fit$b[o], sum=T )
				names(fit$time) = names(a)[o]
				m = match( names(a), names(fit$time))
				fit$time= fit$time[m]
				fit
			}
			fit = NULL
			if ( rev.order ) {
				fit = pred2D( b,a )
				tmp = fit$a
				fit$a= fit$b
				fit$b = tmp
				rm(tmp)
			}else {
				fit = pred2D( a,b )
			}
			
			openPlot( file.path( outpath,"Pseudotime" ) )
			plot( a, b, col= x$usedObj$colorRange[[grouping]][ x$samples[,grouping]] , pch=16)
			o = order( fit$time )
			points( fit$a[o], fit$b[o] , lwd=1.5, col=gplots::bluered(length(a)), pch=16, cex=2 )
			dev.off()
			
			dat = x$dat
			dat@x[which(dat@x == -1)] = 0
			Matrix::drop0(dat)
			
			traj.corGenes = FastWilcoxTest::CorMatrix( dat, fit$time )
			names(traj.corGenes) = rownames(x)
			x$annotation$cor2pseudotime = traj.corGenes 
			traj.corGenes =traj.corGenes[ which( ! is.na(traj.corGenes))]
			
			genes = c( sort(names(traj.corGenes[order(traj.corGenes)[1:n]])), sort(names(traj.corGenes[order(traj.corGenes, decreasing=T)[1:n]]))   ) 
			
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
				unlist(lapply ( (smooth+1):length(x), func, x, smooth ) )
			}
			o = order( fit$time )
			X = roll ( fit$time[o], smooth, 'mean' )
			colV = roll ( as.numeric(x$samples[o,grouping]), smooth , 'mean' )
			colD = roll ( as.numeric(x$samples[o,grouping]), smooth , 'table' )
			col = x$usedObj$colorRange[[grouping]][ round( colV ) ]
			col[which(colD != 1)] = 'black'
			
			## create a summary plot for all genee
			if( ! is.null(summaryPlot)){
				
				openPlot(file.path( outpath, paste(summaryPlot,sep="_",'top') ))
				plot( min(X, na.rm=T), 0, xlim=c(min(X, na.rm=T), max(X, na.rm=T) ), ylim=c(0,1), xlab='pseudotime', ylab='smoothed normalized expression', col='white')
				
				for(  gname in sort(names(traj.corGenes[order(traj.corGenes)[1:n]])) ) {
					Y = roll( dat[gname, o], smooth, 'mean' )
					Y = Y - min(Y)
					Y = Y / max(Y)
					#lines(loess( Y ~ X, se=F , span=0.1),  lwd=1, col= x$usedObj$colorRange[[grouping]][ col ] )
					#browser()
					lines( X, Y, lwd=1, col= 'black' )
					lapply( names(table( col)) , function(n) {
								if ( n != 'black'){
									lines( X[which(col == n)], Y[which(col==n)], lwd=1, col= n )
								}
							})
				}
				dev.off()
				
				openPlot(file.path( outpath, paste(summaryPlot,sep="_",'bottom') ))
				plot(min(X, na.rm=T), 0, xlim=c(min(X, na.rm=T), max(X, na.rm=T) ), ylim=c(0,1) , xlab='pseudotime', ylab='smoothed normalized expression',col='white')
				
				for(  gname in sort(names(traj.corGenes[order(traj.corGenes, decreasing=T)[1:n]])) ) {
					Y = roll( dat[gname, o], smooth, 'mean' )
					Y = Y - min(Y)
					Y = Y / max(Y)
					lines( X, Y, lwd=1, col= 'black' )
					lapply( names(table( col)) , function(n) {
								if ( n != 'black'){
									lines( X[which(col == n)], Y[which(col==n)], lwd=1, col= n )
								}
							})
					#lines(loess( Y ~ X, se=F , span=0.1),  lwd=1, col= x$usedObj$colorRange[[grouping]][ col ] )
				}
				dev.off()
			}
			if ( is.null(summaryPlot)){
				for ( gname in unique(c( plotGenes, genes ) ) ){
					fname = paste(sep="_", gname, "traj.corGenes", smooth,"factor" )
					message(fname)
					
					Y = roll( dat[gname, o], smooth, 'mean' )
					openPlot(file.path( outpath, fname))
					
					plot(min(X, na.rm=T), 0, xlim=c(min(X, na.rm=T), max(X, na.rm=T) ), 
							ylim=c(min(Y), max(Y)) , xlab='pseudotime', ylab='smoothed expression',col='white')
					lines( X, Y, lwd=5, col= 'black' )
					lapply( names(table( col)) , function(n) {
								if ( n != 'black'){
									lines( X[which(col == n)], Y[which(col==n)], lwd=5, col= n )
								}
							})
					#plot( x, y ,  main = gname, xlab='pseudotime', ylab='smoothed expression', type='l', lwd=5 )
					#plot(loess( unlist(Y) ~ unlist(X), se=F , span=0.1),  main = gname, col= col,
					#		xlab='pseudotime', ylab='smoothed expression', type='l', lwd=5 )
					
					dev.off()
				}
			}
			x$samples$Pseudotime = fit$time
			x$samples$PseudotimeX = fit$a
			x$samples$PseudotimeY = fit$b
			
			genes
			
		} )
