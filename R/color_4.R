#' @name colors_4
#' @aliases colors_4,BioData-method
#' @rdname colors_4-methods
#' @docType methods
#' @description Create the colour information for a samples or annotation column
#' @param x the BioData object
#' @param name the name of data column to colour
#' @param colFunc a colour function like default = function(x) {rainbow(x)}
#' @param force forcefully overwrite an existing color entry for this column (default FALSE)
#' @title description of function createRFgrouping_col
#' @export 
if ( ! isGeneric('colors_4') ){ methods::setGeneric('colors_4', ## Name
		function ( x, name,  colFunc = NULL, force=FALSE  ) { 
			standardGeneric('colors_4')
		}
)
}else {
	print ("Onload warn generic function 'colors_4' already defined - no overloading here!")
}

setMethod('colors_4', signature = c ('BioData'),
		definition = function ( x, name, colFunc = NULL, force=FALSE  ) {
			if ( length(name) > 1 ) {
				lapply( name, function(n) { colors_4 (x, n, colFunc,force ) } )
			}else {
			if ( is.null(x$usedObj[['colorRange']])){
				x$usedObj[['colorRange']] = list()
			}
			if ( is.null(colFunc) ){
				colFunc = function(x) {grDevices::rainbow(x)}
			}
			
			mix <- function ( col, l ) {
				if ( l > 10 ){
					if ( l %% 2 == 0 ) {
						col <- c(rbind( col[1:(l/2)], col[(l/2+1):l] ))
					}else {
						col <- c(rbind( col[1:(l/2+1)], c(col[(l/2+1):l],'white') ))[1:l]
					}
				}
				col
			}
			if ( force ) {
				x$usedObj[['colorRange']][[name]] = NULL
			}
			if ( is.null( x$usedObj[['colorRange']][[name]] ) ) {
				if ( !is.na( match(name, colnames(x$samples)))){
					if ( is.na(match(class(x$samples[, name ]),'factor')) ){
						x$samples[, name ] <- factor(x$samples[, name ])
					}
					if ( class(x$samples) == 'factor' ) {
						x$samples = data.frame(x$samples)
					}
					x$usedObj[['colorRange']][[name]] <- mix(
							colFunc( length(levels( x$samples[, name ]))) ,
							length(levels( x$samples[, name ])) 
					)
				}
				else if ( !is.na( match(name, colnames(x$annotation)))){
					if ( is.na(match(class(x$annotation[, name ]),'factor')) ){
						x$annotation[, name ] <- factor(x$annotation[, name ])
					}
					if ( class(x$annotation) == 'factor' ) {
						x$annotation = data.frame(x$annotation)
					}
					x$usedObj[['colorRange']][[name]] <- mix( 
							colFunc( length(levels( x$annotation[, name ]))),
							length(levels( x$annotation[, name ]))
					)
				}
				else {
					stop(paste( "Sorry the column '",name,"' is nether defined in the samples nor in the annotation table!", sep="") )
				}
			}
			}
			invisible(x)
		}
)