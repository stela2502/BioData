#' @name plot_velocyto
#' @aliases plot_velocyto,BioData-method
#' @rdname plot_velocyto-methods
#' @docType methods
#' @description use a velocyto input and plot the velocities on one of the MDS projections
#' @param x The BioData object
#' @param velo the velocyto read reads
#' @param veloCol the column in the BioData samples table containing the velocyto cell names
#' @param group the group to be plotted
#' @param mds_name the mds to plot on (default "Expression PCA")
#' @param ret an optional list gotten from the last plot_velocyto run (default=list())
#' @param ofile an optional outfile for the figure (pdf) default= NULL
#' @param arrow.scale the arrow.scale in the final figure (default 5)
#' @param ... unused at the moment
#' @title description of function plot_velocyto
#' @return a list of velocyto.R results (reusable)
#' @export 
setGeneric('plot_velocyto', ## Name
	function ( x, velo, veloCol, group, mds_name = "Expression PCA", ret=list(), ofile = NULL, arrow.scale=5, ... ) { ## Argumente der generischen Funktion
		standardGeneric('plot_velocyto') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot_velocyto', signature = c ('BioData'),
	definition = function ( x, velo, veloCol, group, mds_name = "Expression PCA", ret=list(), ofile = NULL, arrow.scale=5, ... ) {
	
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
	
	m = match( x$samples[,veloCol], colnames(velo$spliced) )
	matching_BioDataCols =which(is.na(m)==F)
    OKveloCol = m[which(is.na(m)==F)]
	all.equal(x$samples[ matching_BioDataCols,veloCol], colnames(velo$spliced)[OKveloCol] )
	
	velo$spliced   = velo$spliced[,OKveloCol]
	velo$unspliced = velo$unspliced[,OKveloCol]
	velo$ambiguous = velo$ambiguous[,OKveloCol]
		
	colnames(velo$spliced) = colnames(velo$unspliced) = colnames (x$dat )[matching_BioDataCols]
	
	emat <- velo$spliced
	nmat <- velo$unspliced

	# calculate using velocyto.R
	if ( is.null( ret$cell.dist)) {
		print ("calculate cell distances")
		ret$rvel.cd = NULL
		ret$velocity_vals = NULL
		ret$cell.dist <- as.dist(1-velocyto.R::armaCor(t(x$usedObj$pr@scores[matching_BioDataCols,])))
	}
	
	emb = mds(x, mds.type=mds_name, dim=2 )[matching_BioDataCols,]
	#emb <- x$usedObj$MDS_PCA100[[mds_name]][matching_BioDataCols,1:2]
	
	cluster.label=x$samples[matching_BioDataCols,group]
	names(cluster.label) = colnames(x$dat)[matching_BioDataCols]
	
	
	emat <- velocyto.R::filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
	nmat <- velocyto.R::filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
	
	fit.quantile <- 0.02
	if ( is.null(ret$rvel.cd)) {
		print ("calculate cell velocity estimates")
		ret$rvel.cd <- velocyto.R::gene.relative.velocity.estimates(
				emat,nmat,deltaT=1,kCells=20,cell.dist=ret$cell.dist,fit.quantile=fit.quantile)
	}
	
	## now get the colors
	
	#col = x$usedObj$colorRange[[group]][ x$samples[matching_BioDataCols, group] ]
	#names(col) = colnames(x$dat)[matching_BioDataCols]
	
	col = x$usedObj$colorRange[[group]][ x$samples[, group] ]
	names(col) = colnames(x$dat)
	
	## and plot
	print ("create plot")
	## here is where we need to create the plot outfile
	if ( ! is.null(ofile)) {
		png( file=paste(ofile ,'png',sep='.'), width=800, height=800)
	}
	if ( ! is.null(ret$velocity_vals$cc))
	ret$velocity_vals = velocyto.R::show.velocity.on.embedding.cor(
			mds(x, mds.type=mds_name, dim=2 ), ret$rvel.cd, n=300, scale='sqrt', cell.colors=col,
			cex=0.8, arrow.scale=arrow.scale, show.grid.flow=TRUE, min.grid.cell.mass=0.5,
			grid.n=40, arrow.lwd=1, do.par=F, cell.border.alpha = 0.1, return.details=T,
			cc= ret$velocity_vals$cc
	)
	if ( ! is.null(ofile)) {
		dev.off()
	}
	invisible(ret)	
}  )
