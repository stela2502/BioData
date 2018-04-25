#' @name get.genes.cor.to
#' @aliases get.genes.cor.to,BioData-method
#' @rdname get.genes.cor.to-methods
#' @docType methods
#' @description VR function that exports all genes correlating to the input gene name
#' @param x the BioData object
#' @param gname the gene name of interst
#' @param output the VR outpath
#' @title description of function get.genes.cor.to
#' @export 
if ( ! isGeneric('get.genes.cor.to') ){ setGeneric('get.genes.cor.to', ## Name
	function (x,gname,output) { ## Argumente der generischen Funktion
		standardGeneric('get.genes.cor.to') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'get.genes.cor.to' already defined - no overloading here!")
}

setMethod('get.genes.cor.to', signature = c ('BioData'),
	definition = function (x, gname, output) {
	
		if ( length(gname) == ncol(x$data())) {
			print ("You have given me the correlating variables as gname - hope that was intended!")
			goi <- as.numeric(gname)
		}else {
			goi <- as.vector(t(x$data()[gname,]))
		}
	
	
		calc.cor <- function(v,comp){
			cor(v,comp)
		}
	
		cor.values <- apply(x$data(),1,calc.cor,comp=goi)
	
		sort(cor.values)
	
} )

