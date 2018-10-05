
plot_velocyto <- function( x, velo, translation, group, mds_name = "Expression PCA", ret=list(), ofile = NULL, ... ) {
	
	if (!requireNamespace("velocyto.R", quietly = TRUE)) {
		stop("velocyto.R needed for this function to work. Please install it.",
				call. = FALSE)
	}
    # inspired by http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html
	# we have
	# a BioData object with mds and gene selectzion applied.
	# a velocyto object containing the read data
	# and a translation table between velocyto and BioData sample names (matrix with velo and BioData columns)
	# and the group you want to color on.
	
	# What we will do here:
	# get rid of all velocyto samples that are not in our dataset
	
	
	
	ok = match( colnames(velo$spliced), translation$velo )
	ok = which(is.na(ok)==F)
	velo$spliced = velo$spliced[,ok]
	velo$unspliced = velo$sunspliced[,ok]
	velo$ambiguous = velo$ambiguous[,ok]
	
	colnames(velo$spliced) = colnames(velo$unspliced) = 
			translation[match(colnames(velo$spliced), translation[,'velo']  ),'BioData']
	emat <- velo$spliced
	nmat <- velo$unspliced
	
	# calculate using velocyto.R
	if ( is.null( ret$cell.dist)) {
		print ("calculate cell distances")
		ret$rvel.cd = NULL
		ret$velocity_vals = NULL
		ret$cell.dist <- as.dist(1-velocyto.R::armaCor(t(x$usedObj$pr@scores)))
	}
	
	emb <- x$MDS_PCA100[[mds_name]][,1:2]
	
	emat <- velocyto.R::filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
	nmat <- velocyto.R::filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
	
	fit.quantile <- 0.02
	if ( is.null(ret$rvel.cd)) {
		print ("calculate cell velocity estimates")
		ret$rvel.cd <- velocyto.R::gene.relative.velocity.estimates(
				emat,nmat,deltaT=1,kCells=20,cell.dist=ret$cell.dist,fit.quantile=fit.quantile)
	}
	
	## now get the colors
	col = data$usedObj$colorRange[[group]][ data$samples[, group] ]
	names(col) = colnames(x$dat)
	
	## and plot
	print ("create plot")
	## here is where we need to create the plot outfile
	if ( ! is.null(ofile)) {
		pdf( file=paste(ofile ,'pdf',sep='.'), width=8, height=8)
	}
	ret$velocity_vals = velocyto.R::show.velocity.on.embedding.cor(
			emb, ret$rvel.cd, n=300, scale='sqrt', cell.colors=col,
			cex=0.8, arrow.scale=5, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
			grid.n=40, arrow.lwd=1, do.par=F, cell.border.alpha = 0.1, return.details=T
	)
	if ( ! is.null(ofile)) {
		dev.off()
	}
	invisible(ret)	
} 
