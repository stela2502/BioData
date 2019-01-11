#' @name defineHeatmapColors
#' @aliases defineHeatmapColors,BioData-method
#' @rdname defineHeatmapColors-methods
#' @docType methods
#' @description  uses ggplot2 to plot heatmaps
#' @param x the BioData object
#' @param melted the merged data object with the Expression column that should be colored
#' @param colrs and optional colors vector( gray + bluered for data and rainbow for samples)
#' @param lowest the lowest color in the heatmap (gray by default)
#' @title description of function defineHeatmapColors
#' @return a list with the modified merged table and the colors vector
#' @export defineHeatmapColors
if ( ! isGeneric('defineHeatmapColors') ){ setGeneric('defineHeatmapColors', ## Name
		function (x, melted, colrs=NULL, lowest='gray',...) { 
			standardGeneric('defineHeatmapColors')
		}
)
}else {
	print ("Onload warn generic function 'defineHeatmapColors' already defined - no overloading here!")
}

setMethod('defineHeatmapColors', signature = c('tRNAMINT') ,
		definition = function (x, melted,colrs=NULL, lowest='gray', ... ){
			if ( ! is.null(lowest) ) {
				colors= c(
						lowest, 
						gplots::bluered(length(brks) -2  ), ## the expression
				)
			}else {
				colors= c(
						gplots::bluered(length(brks) -1  ), ## the expression
				)
			}
			if ( is.factor( melted$Expression )) {
				## here might be some row grouping going on!
				d <- levels(melted$Expression)[melted$Expression]
				prob.id <- which(is.na(as.numeric(d))==T)
				treat.separate <- unique(d[prob.id])
				n <- as.numeric(d[-prob.id])
				m <- minValueExpr( x)
				brks= c( (m -.1), as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
				brks = unique(brks)
				d[-prob.id]  <- brks [cut( n, breaks= brks)]
				melted$Expression <- factor( d, levels= c(brks, treat.separate ) )
				colors <- c( colors,rainbow( length(treat.separate) ) )
			}
			else {
				n <- as.numeric(melted$Expression )
				m <- minValueExpr( x)
				brks= c( (m-.1), as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
				brks = unique(brks)
				melted$Expression <- factor( brks [cut( n, breaks= brks)] , levels= c(brks) )
			}
			list (melted = melted, colors = colors)
		}
)


setMethod('defineHeatmapColors', signature = c('BioData') ,
		definition = function (x, melted, colrs=NULL, lowest='gray' ){
			if ( is.factor( melted$Expression )) {
				## here might be some row grouping going on!
				d <- levels(melted$Expression)[melted$Expression]
				prob.id <- which(is.na(as.numeric(d))==T)
				treat.separate <- unique(d[prob.id])
				n <- as.numeric(d[-prob.id])
				brks= c( -20.1, as.vector(quantile(n[which(n != -20)],seq(0,1,by=0.1)) ))
				brks = unique(brks)
				d[-prob.id]  <- brks [cut( n, breaks= brks)]
				melted$Expression <- factor( d, levels= c(brks, treat.separate ) )
				colors= c(
						gplots::bluered(length(brks) -1  ), ## the expression
						rainbow( length(treat.separate) ) ## the sample descriptions
				)
			}
			else {
				n <- as.numeric(melted$Expression )
				brks= c(  as.vector(quantile(n,seq(0,1,by=0.1)) ))
				brks = unique(brks)
				melted$Expression <- factor( brks [cut( n, breaks= brks)] , levels= c(brks) )
				
				colors= c(
						lowest, 
						gplots::bluered(length(brks) -2  ) ## the expression
				)
			}
			list (melted = melted, colors = colors)
		}
)

