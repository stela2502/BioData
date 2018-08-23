#' @name as_cellexalvr
#' @aliases as_cellexalvr,cellexalvr-method
#' @rdname as_cellexalvr-methods
#' @docType methods
#' @description Convert a R6::BioData object into a R6::cellexalvr object
#' @param x the BioData object
#' @title description of function as_cellexalvr
#' @export 
setGeneric('as_cellexalvr', ## Name
		function ( x, ... ) { ## Argumente der generischen Funktion
			standardGeneric('as_cellexalvr') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
#' @name as_cellexalvr.BioData
#' @aliases as_cellexalvr.BioData,cellexalvr-method
#' @rdname as_cellexalvr-methods
#' @docType methods
#' @description Convert a R6::BioData object into a R6::cellexalvr object
#' @param x the BioData object
#' @param cellInfo the samples columns containing 0/1 info for the meta.cells matrix
#' @param groups the R calculated groupings to go into the userGroups matrix
#' @param linear linear sample information that is to be stored in the facs matrix
#' @title description of function as_cellexalvr
#' @export 
setMethod('as_cellexalvr', signature = c ('BioData'),
		definition = function ( x, cellInfo, groups=NULL, linear=NULL  ) {
			x <- x$clone()
			print ( "Are you sure, that the rownames of the data object are 'gene symbols'? If not please change that and re-run this function" )
			fixNames <- function( names ) { 
				unlist(lapply( names, 
								function(gname) { paste(unlist(stringr::str_split( gname, '\\s+')), collapse='_') }) )
			} 
			colnames(x$dat) <- make.names( colnames(x$dat) )
			lapply(c ('MDS', 'MDS_PCA100'), function (n) {
					if ( ! is.null(x$usedObj[[n]]) ){
						lapply( names(x$usedObj[[n]]),
								function(v){	rownames(x$usedObj[[n]][[v]]) = colnames(x$dat) }
						)
					}
					else { x$usedObj[[n]] = list()}
				}
			)
			names(x$usedObj[['MDS_PCA100']]) = paste( 'PCA100', names(x$usedObj[['MDS_PCA100']]) )
			
			meta.cell = cellexalvrR::make.cell.meta.from.df( as.matrix(x$samples), cellInfo )
			rownames(meta.cell) <- colnames(x$dat)
			
			userGroups = NULL
			if ( ! is.null(groups) ) {
				newN <- fixNames(groups)
				userGroups = x$samples[,groups]
				colnames(userGroups) = newN
				rownames(userGroups) = colnames(x$dat)
				apply( cbind(groups, newN ),1, 
						function(rnames) { 
							if ( is.null(x$usedObj$colorRange[[ rnames[2] ]])){ 
								if (  is.null(x$usedObj$colorRange[[ rnames[1] ]]) ){ 
									colors_4(x,rnames[2])
								}else {
									x$usedObj$colorRange[[ rnames[2] ]] = x$usedObj$colorRange[[ rnames[1] ]]
								} 
							}  } 
				)
			}
			if ( is.null(userGroups)){
				userGroups = data.frame()
			}
			
			index = NULL
			if ( !is.null(linear) ) {
				index = as.matrix(x$samples[,linear])
				rownames(index) = colnames(x$dat)
				colnames(index) =  fixNames(linear)
			}
			if ( is.null(index) ) {
				index = matrix()
			}
			
			ret = new( 'cellexalvr',
					data = as.matrix(x$dat),
					meta.cell= meta.cell, 
					userGroups= userGroups,
					meta.gene= as.matrix(x$annotation), 
					mds=c ( x$usedObj$MDS, x$usedObj$MDS_PCA100),
					index = index
			)
			names(ret@mds) = unlist(lapply(names(ret@mds), function(n) { str_replace_all(n, "\\s+", '_') } ))
			m <- which(ret@data == -1 )
			if ( length(m) > 0 ) {
				ret@data[m] = 0
			}
			ret@colors = x$usedObj$colorRange
			ret
		}
)
