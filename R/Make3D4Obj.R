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
#' @param size the size of the 3D points default = 3.0
#' @param green single cell normalization looses gene expression values. Display the cells with lost expression in green or black default=FALSE
#' @param useRaw base the projection on the raw data and not the n=100 PCA data (default FALSE)
#' @title description of function Make3D4obj
#' @export 
if ( ! isGeneric('Make3D4obj') ){ methods::setGeneric('Make3D4obj', ## Name
	function ( x, group, mds.type='PCA', cex=0.5, colFunc = function(x) {rainbow(x)}, cut=F, 
			names=F, opath=NULL, main='', genes=F, plotType=1, size=3.0, green=FALSE, useRaw=FALSE ) { 
		standardGeneric('Make3D4obj')
	}
)
}else {
	print ("Onload warn generic function 'Make3D4obj' already defined - no overloading here!")
}


setMethod('Make3D4obj', signature = c ('BioData'),
	definition = function ( x, group, mds.type='PCA', cex=0.5, 
			colFunc = function(x) {rainbow(x)}, cut=F, names=F, opath=NULL, main='', genes=F , plotType=1, size=3.0 , green=FALSE, useRaw=FALSE ) {

		MDS_NAME = 'MDS_PCA100'
		My.legend3d <- function (main= '', ...) {
			if ( ! exists ( 'main')) {
				main = ''
			}
			rgl::bgplot3d( {
				graphics::par( mar =c(1,1,1,1),bg='#4C4C4C')
				graphics::plot(0, 0, type = "n", xlim = 0:1, ylim = 0:1, xaxs = "i",
						yaxs = "i", axes = FALSE, bty = "n", col='#4C4C4C', main=main, col.main =  "white" )
				graphics::legend(...)
			} )
		}
		check_and_replace <- function( name, list) {
			if ( length(grep(name, names(list) )) ==1 ){
				name = names(list)[grep(name, names(list) )]
				print ( paste( "name changed to",name))
			}
			name
		}
		if ( genes ) {
			x = transpose(x$clone())
			return (Make3D4obj(x,group=group, mds.type=mds.type,  cex=cex, colFunc=colFunc, cut= cut, names=names,
							opath=opath, main=main, genes=F, plotType=plotType, size=size, green=green, useRaw=useRaw));
		}
		MDS_NAME = 'MDS_PCA100'
		if ( useRaw ){
			MDS_NAME = 'MDS'
		}
		mds.type = check_and_replace( mds.type, x$usedObj[[MDS_NAME]] )

		if ( is.null (x$usedObj[[MDS_NAME]][[mds.type]] )){
			mds(x, mds.type=mds.type, useRaw = useRaw)
			mds.type = check_and_replace( mds.type, x$usedObj$MDS )
		}
		
		title = paste( group, mds.type )
		
        if ( cut ) {
                ## this is a gene expression value!
                n <- as.numeric(x$data()[group,] )
                m <- min( n )
                
				col = NULL
				COLS = NULL
				if ( m == -1 | m == -21){
					brks= c( (m-.1),m+1-.1 , m+1 ,as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
              	  	brks = unique(as.numeric(sprintf("%2.6e", brks)))
					d  <- factor(brks [cut( n, breaks= brks)], levels=brks)
					if ( green) {
						COLS = c('#006D2C', 'black', gplots::bluered(length(brks) -1  ))
					}else{
						COLS = c('black', 'black', gplots::bluered(length(brks) -1  ))
					}
					
					col = COLS[d]
				}else {
					brks= c( (m-.1),m,as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
					brks = unique(as.numeric(sprintf("%2.6e", brks)))
					d  <- factor(brks [cut( n, breaks= brks)], levels=brks)
					COLS = c('black', gplots::bluered(length(brks) -1  ))
					col = COLS[d]
				}      
                
        }else {
			colors_4(x, group, colFunc=colFunc )
			col <- x$usedObj$colorRange[[group]][x$samples[,group]]
		}
		#18 105 762 810
		rgl::rgl.open()
		rgl::par3d(windowRect = c(18,105, 762, 810))
		Sys.sleep(1)
		#bg3d(color='#4C4C4C') 
		if ( plotType == 1) {
        if ( cut ) {
                ## plot points!
				#browser()
				My.legend3d ("topright", legend = paste( brks ), pch=16, col= COLS, cex=1,inset =c(0.02),
						main = title
				)
				mds.type = check_and_replace( mds.type, x$usedObj[[MDS_NAME]] )
				if ( is.null(x$usedObj[[MDS_NAME]][[mds.type]])){
					MDS_NAME = stringr::str_replace(MDS_NAME, '_PCA100', '' )
				}
                rgl::rgl.points( x$usedObj[[MDS_NAME]][[mds.type]], col=col, size=size )

        }
        else {
                if ( names ) {
						My.legend3d ("topright", legend = paste( unique(levels(x$samples[,group]))  ), pch = 16, 
								col = x$usedObj$colorRange[[group]], cex=1, inset=c(0.02),
								main = title
						)
						mds.type = check_and_replace( mds.type, x$usedObj[[MDS_NAME]] )
						if ( is.null(x$usedObj[[MDS_NAME]][[mds.type]])){
							MDS_NAME = stringr::str_replace(MDS_NAME, '_PCA100', '' )
						}
                        rgl::rgl.texts( x$usedObj[[MDS_NAME]][[mds.type]], col=col, text= as.character(x$samples[,group]), cex=cex )
                }
                else {
                        My.legend3d ("topright", legend = paste( unique(levels(x$samples[,group]))  ), 
								pch = 16, col = x$usedObj$colorRange[[group]], cex=1, inset=c(0.02),
								main = title
						)
						mds.type = check_and_replace( mds.type, x$usedObj[[MDS_NAME]] )
						if ( is.null(x$usedObj[[MDS_NAME]][[mds.type]])){
							MDS_NAME = stringr::str_replace(MDS_NAME, '_PCA100', '' )
						}
						rgl::rgl.points( x$usedObj[[MDS_NAME]][[mds.type]], col=col, size=size )
						
                }
        }
		}else if ( plotType == 2) {
			rgl::bg3d("white")
			if ( names) {
				rgl::rgl.texts( x$usedObj[[MDS_NAME]][[mds.type]], col=col, text= as.character(x$samples[,group]), cex=cex )
			}else {
				rgl::rgl.points( x$usedObj[[MDS_NAME]][[mds.type]], col=col, size=size )
			}
			rgl::grid3d(c("x", "y", "z"))
			rgl::axis3d(c("x+"),col="black",xlab="Component 1")
			rgl::axis3d(c("y+"),col="black")
			rgl::axis3d(c("z+"),col="black")
			
		}
} )


