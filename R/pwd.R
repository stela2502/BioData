#' @name pwd
#' @aliases pwd
#' @rdname pwd-methods
#' @docType methods
#' @description  uses the linux pwd command to determin the working directory 
#' @return A string containing the working directory 
#' @title description of function pwd
#' @export 
setGeneric('pwd', ## Name
	function ( a ) { ## Argumente der generischen Funktion
		standardGeneric('pwd') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('pwd', signature = c () ,
	definition = function ( a ) {
		rm(a)
	system( 'pwd > __pwd' )
	t <- read.delim( file = '__pwd', header=F)
	t <- as.vector(t[1,1])
	t <- paste(t,"/",sep='')
	unlink( '__pwd')
	t
})

