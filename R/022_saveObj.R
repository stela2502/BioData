#' BioData has moved from, being my main analysis pipeline to mere herlper status.
#' 
#' I do not wnat to remove anything here, but instead implement a way to handle Seurat objects, too.
#' 
#' @name saveObj
#' @aliases saveObj,ANY-method
#' @rdname saveObj-methods
#' @docType methods
#' @description This function saves the object either as analysis.RData or norm_data.RData if the analysi.RData has not been produced before
#' @param data the BioData, Seurat or something else object
#' @param file the outfile default=paste(str_replace_all(data$name, "\\s+", "_" ), '.RData',sep='')
#'        the file parameter has to be set for all but BioData objects!
#' @title description of function saveObj
#' @export 
if ( ! isGeneric('saveObj') ){ methods::setGeneric('saveObj', ## Name
		function ( data, file=NULL ){	
			standardGeneric('saveObj')
		}
)
}else {
	print ("Onload warn generic function 'saveObj' already defined - no overloading here!")
}

setMethod('saveObj', signature = c ('ANY'),
		definition = function ( data, file=NULL ){
			if ( is.null(file)){
				file=paste(stringr::str_replace_all(data$name, "\\s+", "_" ), '.RData',sep='')
			}
			path = NULL
			if ( length (grep('BioData', class(data))) != 0 ){
				path = data$outpath
				data$dat = Matrix::drop0(data$dat)
				if ( !is.null(data$raw)) {
					data$raw = Matrix::drop0(data$raw)
				}
				if ( !is.null(data$zscored)) {
					data$zscored = Matrix::drop0(data$zscored)
				}
			}else {
				if ( is.null(file)) {
					stop("None BioData objects need to have a file to save to")
				}
				path = getwd()
			}
			print ( paste('data exported to', file.path(path,file) ) )
			print ( timestamp(quiet=T) )
			gc()
			save(data , file=file.path(path, file) )
			print ({ts=tools::md5sum( file.path(path, file))} )
			
		}
)


