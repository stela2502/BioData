#' @name loadObj
#' @aliases loadObj,BioData-method
#' @rdname loadObj-methods
#' @docType methods
#' @description This function loads the a BioData file and returns the data R object int them.
#'  The function is mainly to keep the scripts simpler. There is no magic involved.
#' @param file the infile
#' @title description of function loadObj
#' @export 
setGeneric('loadObj', ## Name
		function ( file=NULL ){	## Argumente der generischen Funktion
			standardGeneric('loadObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('loadObj', signature = c ('character'),
		definition = function ( file=NULL ){
			if ( is.null(file)){
				stop( "Sorry I need a 'file' to load from" )
			}
			data <- NULL # will be read by the next line
			load( file )
			data
		}
)


