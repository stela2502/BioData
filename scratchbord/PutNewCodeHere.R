pseudotimeTest <- function ( x, a, b, grouping, outpath=NULL,  n=100, plotGenes=NULL ) {

	if ( is.null( outpath )) {
		outpath = x$outpath
	}
	fit = lm( b~a )
	
	png ( file="Pseudotime.png", width=800, height=800)
	plot( a ,  b, col= x$usedObj$colorRange[[grouping]][ x$samples[,grouping]] )
	lines( a,  fit$fitted.values, lwd=1.5, col='red' )
	dev.off()
	dat = x$dat
	dat@x[which(dat@x == -1)] = 0
	traj.corGenes = FastWilcoxTest::CorMatrix( dat, fit$fitted.values )
	names(traj.corGenes) = rownames(x)
	traj.corGenes =traj.corGenes[ which( ! is.na(traj.corGenes))]
	
	genes = c( sort(names(traj.corGenes[order(traj.corGenes)[1:n]])), sort(names(traj.corGenes[order(traj.corGenes, decreasing=T)[1:n]]))   ) 
	
	for ( gname in unique(c( plotGenes, genes ) ) ){
		fname = paste(sep="_", gname, "traj.corGenes", smooth,"factor.png" )
		message(fname)
		d = as.data.table(list( fit$fitted.values, dat[gname, ] ))
		d[which(d[,'V2'] < 0), 'V2'] = 0
		png( file.path( 'Tail_smoothed_expression', fname), width=800, height=800)
	
		x = frollmean(d[, V1], smooth)
		y = frollmean(d[, V2], smooth)
		plot( x, y ,  main = gname, xlab='pseudotime', ylab='smoothed expression', type='l', lwd=5 )
#       plot(loess( y ~ x, se=F , span=0.1),  main = gname, xlab='pseudotime', ylab='smoothed expression', type='l', lwd=5 )
		dev.off()
	}
	
}

