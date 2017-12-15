
#' @name complexHeatmap
#' @aliases complexHeatmap,BioData-method
#' @rdname complexHeatmap-methods
#' @docType methods
#' @description plot the PCR heatmap using the heatmap.3 function included in this package 
#' @param x the BioData object
#' @param ofile the outfile to create in the x$outpath folder
#' @param colGroups columns in the samples table to use to order the data (first == order)
#' @param rowGroups rows in the annotation table to use to color the heatmap rows (first == order)
#' @param colColors a named list of column color vectors
#' @param rowColors a named list of row color vectors
#' @param pdf export as pdf (default = FALSE)
#' @param subpath the subpath for the plots (default = '')
#' @param heapmapCols the color function to calculate the heatmap colours ( default function (x) { c("darkgrey",bluered(x)) } )
#' @param brks how many breaks should the expression value color key have (default=10)
#' @title description of function complexHeatmap
#' @export 
setGeneric('complexHeatmap', ## Name
		function ( x,  ofile=NULL, colGroups=NULL, rowGroups=NULL, colColors=NULL, rowColors=NULL, pdf=FALSE, subpath='', 
				main = '',  heapmapCols= function(x){ c("darkgrey",bluered(x))}, brks=10, X11type= 'cairo' ) { ## Argumente der generischen Funktion
			standardGeneric('complexHeatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('complexHeatmap', signature = c ('BioData'),
		definition = function ( x,  ofile=NULL, colGroups=NULL, rowGroups=NULL, colColors=NULL, rowColors=NULL, pdf=FALSE,
				subpath='', main = '' ,  heapmapCols= function(x){ c("darkgrey",bluered(x))}, brks=10, X11type= 'cairo' ) {
			
			Rowv = FALSE
			Colv = FALSE
			dendrogram = 'both'
			ColSideColors <- NULL
			RowSideColors <- NULL
			ColSideColorsSize <- 1
			RowSideColorsSize <- 1
			x <- x$clone()
			if ( is.null(colColors) ){
				colColors <- x$usedObj[['colorRange']]
			}
			if ( is.null(rowColors) ){
				rowColors <- x$usedObj[['colorRange']]
			}
			if ( ! is.null(colGroups) ) {
				ColSideColorsSize <- length(colGroups)
				x <- reorder.samples(x, colGroups[1] )
				for ( i in colGroups ){
					if ( is.na(match( i, names(colColors))) ){
						x <- colors_4( x, i )
						colGroups[[i]] <- x$usedObj[['colorRange']][[i]]
						#stop( paste( "No colours for the grouping", i, "in the colour objects:", paste(names(colColors), collapse= ' , ') ) )
					}
					ColSideColors <- cbind(ColSideColors, colColors[[ match( i, names(colColors)) ]][x$samples[, i]] )
				}
				colnames(ColSideColors) = colGroups
				#ColSideColors <- matrix( ColSideColors, ncol= ColSideColorsSize)
				Colv = FALSE
				if ( !is.null(rowGroups)){
					dendrogram = 'none'
				}else{
					dendrogram= 'none'
				}
			}else {
				## probably calculate the clustering??
			}
			if ( ! is.null(rowGroups) ) {
				RowSideColorsSize <- length(rowGroups)
				x <- reorder.genes(x, rowGroups[1] )
				for ( i in rowGroups ){
					if ( is.na(match( i, names(rowColors))) ){
						x <- colors_4( x, i )
						rowColors[[i]] <- x$usedObj[['colorRange']][[i]]
						#stop( paste( "No colours for the grouping", i, "in the colour objects:", paste(names(colColors), collapse= ' , ') ) )
					}
					RowSideColors <- rbind( RowSideColors,rowColors[[ match( i, names(rowColors)) ]][x$annotation[, i]] )
				}
				rownames(RowSideColors) = rowGroups
				Rowv = FALSE
				#RowSideColors <- matrix( RowSideColors, nrow= RowSideColorsSize)
				if ( !is.null(colGroups)){
					dendrogram = 'none'
				}else{
					dendrogram= 'none'
				}
			}else {
				## probably calculate the clustering??
			}
			data <- as.matrix(x$data())
			
			brks <- unique(as.vector(c(min(data),quantile(data[which(data!= min(data))],seq(0,1,by=1/brks)),max(data))))
			if ( ! is.null(ofile)){
				if ( pdf ) {
					width= ceiling(nrow(x$samples)/300) * 10
					height = ceiling( nrow(x$annotation) / 100 ) * 10
					if ( height < 8){
						height = 8
					}
					pdf( file=paste(file.path(x$outpath,ofile),'pdf',sep='.'), width=10, height=height)
				}else{
					width= ceiling(nrow(x$samples)/300) * 1600
					height = ceiling( nrow(x$annotation) / 100 ) *800
					png( file=paste(file.path(x$outpath,ofile),'png',sep='.'), width=1600, height=800, type=X11type)
				}
				for ( v in colGroups ) {
					plotLegend(x, file=paste(ofile, 'col'), colname=v, pdf=pdf, col=colColors[[v]], X11type=X11type )
				}
				for ( v in rowGroups ) {
					plotLegend(x, file=paste(ofile, 'row'), colname=v, pdf=pdf, col=rowColors[[v]], X11type=X11type )
				}
			}
			heatmap.3(
					data, breaks=brks,col=heapmapCols(length(brks)-2), Rowv= is.null(RowSideColors), Colv = is.null(ColSideColors),  key=F, symkey=FALSE,
					trace='none', 
					ColSideColors=ColSideColors,ColSideColorsSize=ColSideColorsSize, 
					RowSideColors=RowSideColors,RowSideColorsSize=RowSideColorsSize, 
					cexRow=0.6,cexCol=0.7,main=main, dendrogram=dendrogram, labCol = "", lwid=c(0.5,4), lhei=c(1,4)
			)
			
			if ( ! is.null(ofile)){
				dev.off()
				fn <- paste(file.path(x$outpath,x$name),'_legend_values.pdf',sep='.')
				if ( ! file.exists(fn) ){
				pdf( file=fn, width=8, height=4)
				Z <- as.matrix(1:(length(brks)-2))
				image(Z, col=heapmapCols(length(brks)-2),axes = FALSE, main='color key')
				axis( 1, at=c(0,0.1,1), labels=c('NA','low','high'))
				dev.off()
				}
			}
			
		}
)
