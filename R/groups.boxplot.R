#' @name groups.boxplot
#' @aliases groups.boxplot,BioData-method
#' @rdname groups.boxplot-methods
#' @docType methods
#' @description  This function can be used to get an overview of the different gene level groups in
#' @description  an BioData It uses the linux montage command to merge all different boxplots
#' @description  into one figure
#' @param x the ExpresionSet object
#' @param groupCol the column in the samples table that contains the sample grouping information
#' @param groupRow the column in the annotation table hat contains the gene grouping information
#' @param svg create the single boxplots as svg or (default png) files
#' @param fname the filename extension files are created by paste(x$outpath,fname,<GroupID>,"_boxplot")
#' @param Collapse if set will lead to a collapse of the sample groups into one value per gene. Supports all \code{\link{collaps}} by options.
#' @title description of function groups.boxplot
#' @export 
if ( ! isGeneric('groups.boxplot') ){ setGeneric('groups.boxplot', ## Name
	function ( x, groupCol='GroupName', clusters, svg=F, fname='group_', width=800, height=800,mar=NULL, Collapse=NULL, ...) { ## Argumente der generischen Funktion
		standardGeneric('groups.boxplot') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'groups.boxplot' already defined - no overloading here!")
}

setMethod('groups.boxplot', signature = c ( 'BioData') ,
	definition = function ( x, groupCol, groupRow, svg=F, fname='group_', width=800, height=800,mar=NULL, Collapse=NULL, ...) {
		clusters <- as.numeric( x$annotation[,groupRow])
		maxG <- max( clusters )
	ret <- list()
	r=1
	gnames <- names(table(x$samples[,groupCol]))
	if ( ! is.null( Collapse) ){
		x <- collaps( x, by = Collapse, groupCol = groupCol )
		fname = paste( fname, Collapse,"_",sep='')
	}
	x <- z.score(x)
	fnames <- vector('character', maxG )
	ma = -100
	mi = +100
	for ( i in 1:maxG ){
		if ( svg ) {
			fnames[i] = file.path(x$outpath,paste(fname,i,"_boxplot.C.svg",sep=''))
			devSVG ( file= fnames[i], width=width/130, height=height/130 )
		}else{
			fnames[i] =paste(x$outpath,fname,i,"_boxplot.png",sep='')
			png( file=fnames[i],width=width,height=height )
		}
		
		robj <- reduceTo( x,what='row',to= rownames(x$data())[which(clusters==i)], name=paste("group_",i,sep=''), copy=T )
		ret[[r]] <- rownames(robj$data)
		r = r+1
		names(ret)[length(ret)] = robj$name
		a= 1;
		d <- list()
		for ( n in as.vector(gnames) ){
			d[[a]] = as.vector(unlist(robj$data()[,which(robj$samples[,groupCol] == n)]))
			a = a+1
		}
		names(d) <- gnames
		A <- NULL
		if ( ! is.null(mar) ){
			A <- boxplot(d, main=paste('Cluster', i, ' - ',nrow(robj$data())," genes", sep='' ),outline=FALSE, par(mar=mar), ...  )
		}else {
			A <- boxplot(d, main=paste('Cluster', i, ' - ',nrow(robj$data())," genes", sep='' ),outline=FALSE, ...  )
		}
		mi <- min(c(mi,A$stats[1,]))
		ma <- max(c(ma, A$stats[5,]))
		dev.off()
	}
	print ( paste(  "min",mi, 'max', ma) )
	#print (paste('montage', paste(fnames, collapse= " "), "-geometry", paste("'", width, "x", height,"+0+0'",sep=''),paste(x$outpath,fname,"montage.png",sep=''), sep=' ' ))
	try( file.remove(  file.path(x$outpath,paste(fname,"montage.png",sep='')) ) , silent=T )
	system ( paste('montage', fnames, "-geometry", paste("'", width, "x", height,"+0+0'",sep=''),file.path(x$outpath,paste(fname,"montage.png",sep=''))," 2>/dev/null", collapse=' ' ) )
	ret
})

