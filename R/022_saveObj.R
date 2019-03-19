#' @name saveObj
#' @aliases saveObj,BioData-method
#' @rdname saveObj-methods
#' @docType methods
#' @description This function saves the object either as analysis.RData or norm_data.RData if the analysi.RData has not been produced before
#' @param data the BioData object
#' @param file the outfile default=paste(str_replace_all(data$name, "\\s+", "_" ), '.RData',sep='')
#' @title description of function saveObj
#' @export 
if ( ! isGeneric('saveObj') ){ setGeneric('saveObj', ## Name
		function ( data, file=NULL ){	
			standardGeneric('saveObj')
		}
)
}else {
	print ("Onload warn generic function 'saveObj' already defined - no overloading here!")
}

setMethod('saveObj', signature = c ('BioData'),
		definition = function ( data, file=NULL ){
			if ( is.null(file)){
				file=paste(stringr::str_replace_all(data$name, "\\s+", "_" ), '.RData',sep='')
			}
			print ( paste('data exported to', file.path(data$outpath,file) ) )
			data$dat = drop0(data$dat)
			if ( !is.null(data$raw)) {
				data$raw = drop0(data$raw)
			}
			if ( !is.null(data$zscored)) {
				data$zscored = drop0(data$zscored)
			}
			save(data , file=file.path(data$outpath, file) )
			
		}
)


