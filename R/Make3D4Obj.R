#' @name Make3D4obj
#' @aliases Make3D4obj,BioData-method
#' @rdname Make3D4obj-methods
#' @docType methods
#' @description This function  created a MDS 3D plot with an inbuilt background describing the data colors 
#' @param x  the BioData object
#' @param group the grouping value (either colname from the samples table or a rowname from the data table
#' @param mds.type the mds type to use (names from names(x$usedObj$MDS)) default = PCA
#' @param cex  define the size of the strings default=0.5
#' @param colFunc the color function if the color is not already defined using colors_4()  default= function(x){rainbow(x)}
#' @param cut this has to be true for genes, as the samples are then binned into 10 expression groups each containing the same number of samples default=F
#' @param names if true not dots, but the name in the samples table is plotted in 3D default=F
#' @param opath create a webgl representation of the plot in a html page in this path (default = NULL) 
#' @param main the title of the plot (default ='')
#' @param genes use gene level MDS data not sample level (default = FALSE)
#' @param plotType choose one [1,2] and check whether you like it ;-) default=1 
#' @title description of function Make3D4obj
#' @export 
if ( ! isGeneric('Make3D4obj') ){ setGeneric('Make3D4obj', ## Name
	function ( x, group, mds.type='PCA', cex=0.5, colFunc = function(x) {rainbow(x)}, cut=F, names=F, opath=NULL, main='', genes=F, plotType=1 ) { ## Argumente der generischen Funktion
		standardGeneric('Make3D4obj') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'Make3D4obj' already defined - no overloading here!")
}


setMethod('Make3D4obj', signature = c ('BioData'),
	definition = function ( x, group, mds.type='Expression PCA', cex=0.5, colFunc = function(x) {rainbow(x)}, cut=F, names=F, opath=NULL, main='', genes=F , plotType=1 ) {

		My.legend3d <- function (...) {
			if ( ! exists ( 'main')) {
				main = ''
			}
			bgplot3d( {
				par( mar =c(1,1,1,1),bg='#4C4C4C')
				plot(0, 0, type = "n", xlim = 0:1, ylim = 0:1, xaxs = "i",
						yaxs = "i", axes = FALSE, bty = "n", col='#4C4C4C', main=main)
				legend(...)
			} )
		}
		check_and_replace <- function( name, list) {
			if ( length(grep(name, names(list) )) ==1 ){
				name = names(list)[grep(name, names(list) )]
			}
			name
		}
		if ( genes ) {
			mds.type = check_and_replace( mds.type, x$usedObj$MDSgenes )
			if ( is.null (x$usedObj$MDSgenes[[mds.type]] )){
				x <- mds(x, mds.type=mds.type, genes=T)
			}
			mds.type = check_and_replace( mds.type, x$usedObj$MDSgenes )
		}else {
			mds.type = check_and_replace( mds.type, x$usedObj$MDS )
			if ( is.null (x$usedObj$MDS[[mds.type]] )){		
				x <- mds(x, mds.type=mds.type)
			}
			mds.type = check_and_replace( mds.type, x$usedObj$MDS )
		}
		
		if ( genes ) {
			x <- x$clone()
			t <- transpose(x)
		}
        if ( cut ) {
                ## this is a gene expression value!
                n <- as.numeric(x$data()[group,] )
                m <- min( n )
                brks= c( (m-.1),m ,as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
                brks = unique(as.numeric(sprintf("%2.6e", brks)))
                d  <- factor(brks [cut( n, breaks= brks)], levels=brks)
                col = c('black', bluered(length(brks) -1  ))[d]
        }else {
			colors_4(x, group, colFunc=colFunc )
			col <- x$usedObj$colorRange[[group]][x$samples[,group]]
		}
		#18 105 762 810
		rgl.open()
		par3d(windowRect = c(18,105, 762, 810))
		Sys.sleep(1)
		#bg3d(color='#4C4C4C') 
		if ( plotType == 1) {
        if ( cut ) {
                ## plot points!
				My.legend3d ("topright", legend = paste( brks ), pch=16, col= c('black', bluered(length(brks) -1  )), cex=1,inset =c(0.02))
                rgl.points( x$usedObj$MDS[[mds.type]], col=col )

        }
        else {
                if ( names ) {
						My.legend3d ("topright", legend = paste( unique(as.character(x$samples[,group]))  ), pch = 16, col = unique(col), cex=1, inset=c(0.02))
                        rgl.texts( x$usedObj$MDS[[mds.type]], col=col, text= as.character(x$samples[,group]), cex=cex )
                }
                else {
                        My.legend3d ("topright", legend = paste( unique(as.character(x$samples[,group]))  ), pch = 16, col = unique(col), cex=1, inset=c(0.02))
						rgl.points( x$usedObj$MDS[[mds.type]], col=col )
						
                }
        }
		}else if ( plotType == 2) {
			bg3d("white")
			if ( names) {
				rgl.texts( x$usedObj$MDS[[mds.type]], col=col, text= as.character(x$samples[,group]), cex=cex )
			}else {
				rgl.points( x$usedObj$MDS[[mds.type]], col=col )
			}
			grid3d(c("x", "y", "z"))
			axis3d(c("x+"),col="black",xlab="Component 1")
			axis3d(c("y+"),col="black")
			axis3d(c("z+"),col="black")
			
		}
} )


