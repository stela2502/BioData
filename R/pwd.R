#' @name pwd
#' @aliases pwd
#' @rdname pwd-methods
#' @docType methods
#' @description  uses the linux pwd command to determin the working directory 
#' @return A string containing the working directory 
#' @title description of function pwd
#' @export 
setGeneric('pwd', ## Name
	function ( a ) { 
		standardGeneric('pwd')
	}
)


setMethod('pwd', signature = c () ,
	definition = function ( a ) {
		rm(a)
	getwd()
})

