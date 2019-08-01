
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
#' @param main the picture title (default '')
#' @param heapmapCols the color function to calculate the heatmap colours ( default function (x) { c("darkgrey",bluered(x)) } )
#' @param brks how many breaks should the expression value color key have (default=10)
#' @param X11type sometimes needed for compatibility (default = 'cairo')
#' @param green if in SingleCell mode normalization losses get a -1 and can be displayed as green or black in the default coloring (default green=F => black)
#' @param noBreaks show linear data colors instead of bined ones; impossible if green=TRUE (default =FALSE)
#' @title description of function complexHeatmap
#' @export 
setGeneric('complexHeatmap', ## Name
		function ( x,  ofile=NULL, colGroups=NULL, rowGroups=NULL, colColors=NULL, rowColors=NULL, pdf=FALSE, subpath='', 
				main = '',  heapmapCols= function(x){ c("darkgrey", gplots::bluered(x))}, brks=10, X11type= 'cairo', green = F, noBreaks = F ) { 
			standardGeneric('complexHeatmap')
		}
)

setMethod('complexHeatmap', signature = c ('BioData'),
		definition = function ( x,  ofile=NULL, colGroups=NULL, rowGroups=NULL, colColors=NULL, rowColors=NULL, pdf=FALSE,
				subpath='', main = '' ,  heapmapCols=NULL, brks=10, X11type= 'cairo', green = F, noBreaks = F ) {
			
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
				dendrogram='col'
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
				if ( dendrogram== 'col') {
					dendrogram= 'both'
				}else {
					dendrogram= 'row'
				}
			}
			data <- as.matrix(x$data())
			m <- min(data)
			if ( m == -1) {
				## new z.score function to keep lost info
				brks <- unique(as.vector(c(m, m+1, stats::quantile(data[which(data > m +1 )],seq(0,1,by=1/brks)),max(data))))
				brks[1:2] = brks[1:2] - 1e-4
				if ( is.null(heapmapCols)){
					if ( green ) {
						heapmapCols = function(x){ c("#006D2C","black", gplots::bluered(x-1))}
					}else {
						heapmapCols = function(x){ c("black","black", gplots::bluered(x-1))}
					}
				}
			}else {
				brks <- unique(as.vector(c(m, stats::quantile(data[which(data!= m)],seq(0,1,by=1/brks)),max(data))))
				if ( is.null(heapmapCols)){heapmapCols = function(x){ c("black", gplots::bluered(x))}}
			}
			if ( ! is.null(ofile)){
				## here I need more magic
				x$name = paste( collapse='_',unlist(strsplit( x$name, '\\s+', perl=T)))
				if ( length(grep(.Platform$file.sep, ofile)) == 0 ) {
					ofile <- file.path(x$outpath,ofile) ## you can also put it specifily somewhere else.
				}else {
					x$outpath <- dirname(ofile) ## in case output should go somewhere else
				}
				if ( length(grep( x$name, ofile) ) == 0 ) {
					ofile = file.path( dirname(ofile), paste(x$name, basename(ofile), sep="_") )
				}
				if ( pdf ) {
					width= ceiling(nrow(x$samples)/300) * 10
					height = ceiling( nrow(x$annotation) / 100 ) * 10
					if ( height < 8){
						height = 8
					}
					grDevices::pdf( file=paste(ofile ,'pdf',sep='.'), width=10, height=height)
				}else{
					width= ceiling(nrow(x$samples)/300) * 1600
					height = ceiling( nrow(x$annotation) / 100 ) *800
					grDevices::png( file=paste(ofile,'png',sep='.'), width=1600, height=800, type=X11type)
				}
				for ( v in colGroups ) {
					plotLegend(x, file=paste(x$name, 'col'), colname=v, pdf=pdf, col=colColors[[v]], X11type=X11type )
				}
				for ( v in rowGroups ) {
					plotLegend(x, file=paste(x$name, 'row'), colname=v, pdf=pdf, col=rowColors[[v]], X11type=X11type )
				}
			}
			if ( length(brks) < 3 ) {
				## not good!
				print ("Highly invariant data - your're shure that heatmap is what you want?" )
				if ( is.null(data$zscored) ){
					brks= c( -1, 0, 1, 2, max(x$data()) )
				}
				else {
					brks= c( -3,-2,-1, 0, 1, 2, 3 )
				}
			}
			if (noBreaks) {
				OK = which(data > 0 )
				#browser()
				data[OK] = data[OK] - min(data[OK])+1e-6
				data[which(data>3)] = 3
				data[which(data < 0 )] = 0
				
				heatmap.3(
						data, col=c('black', gplots::bluered(8) ), Rowv= is.null(RowSideColors), Colv = is.null(ColSideColors),  key=T, symkey=FALSE,
						trace='none', 
						ColSideColors=ColSideColors,ColSideColorsSize=ColSideColorsSize, 
						RowSideColors=RowSideColors,RowSideColorsSize=RowSideColorsSize, 
						cexRow=0.6,cexCol=0.7,main=main, dendrogram=dendrogram, labCol = "", 
						lwid=c(0.5,4), lhei=c(1,4)
				)
				if ( ! is.null(ofile)){
					grDevices::dev.off()
				}
			}
			else {
				heatmap.3(
						data, breaks=brks,col=heapmapCols(length(brks)-2), Rowv= is.null(RowSideColors), Colv = is.null(ColSideColors),  key=F, symkey=FALSE,
						trace='none', 
						ColSideColors=ColSideColors,ColSideColorsSize=ColSideColorsSize, 
						RowSideColors=RowSideColors,RowSideColorsSize=RowSideColorsSize, 
						cexRow=0.6,cexCol=0.7,main=main, dendrogram=dendrogram, labCol = "", 
						lwid=c(0.5,4), lhei=c(1,4)
				)
			
				if ( ! is.null(ofile)){
					grDevices::dev.off()
					fn <- paste(file.path(x$outpath,x$name),'_legend_values.pdf',sep='.')
					if ( ! file.exists(fn) ){
						grDevices::pdf( file=fn, width=8, height=4)
						Z <- as.matrix(1:(length(brks)-2))
						graphics::image(Z, col=heapmapCols(length(brks)-2),axes = FALSE, main='color key')
						if ( min(x$data()) == -1) {
						graphics::axis( 1, at=c(0,0.1,0.2,1), labels=c('lost','NA','low','high') )
						}else {
							graphics::axis( 1, at=c(0,0.1,1), labels=c('NA','low','high') )
						}
						grDevices::dev.off()
					}
				}
			}
			invisible(x)
		}
)
