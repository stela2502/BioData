#' @name show
#' @aliases show,BioData-method
#' @rdname show-methods
#' @docType methods
#' @description  show the BioData
#' @param x the BioData object
#' @return nothing
#' @title description of function show
#' @export 
setMethod('show', signature(object=c('BioData')) ,
	definition = function (object) {
	cat (paste("An object of class", class(object)),"\n" )
	cat("named ",object$name,"\n")
	cat (paste( 'with',nrow(object$dat),'genes and', ncol(object$dat),' samples.'),"\n")
	cat (paste("Annotation datasets (",paste(dim(object$annotation),collapse=','),"): '",paste( colnames(object$annotation ), collapse="', '"),"'  ",sep='' ),"\n")
	cat (paste("Sample annotation (",paste(dim(object$samples),collapse=','),"): '",paste( colnames(object$samples ), collapse="', '"),"'  ",sep='' ),"\n")
	if ( length(names(object$stats)) > 0 ){
		cat ( "P values were calculated for ", length(names(object$stats)) -1, " condition(s)\n")
	}
})


