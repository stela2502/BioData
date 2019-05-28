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
if ( ! isGeneric('get.genes.cor.to') ){ methods::setGeneric('get.genes.cor.to', ## Name
	function (x,gname,output) { 
		standardGeneric('get.genes.cor.to')
	}
)
}else {
	print ("Onload warn generic function 'get.genes.cor.to' already defined - no overloading here!")
}

setMethod('get.genes.cor.to', signature = c ('BioData'),
	definition = function (x, gname, output) {
	
		
		this = x$data()
		bad = which(this@x == -1 )
		if ( length(bad > 0 ) )
			this@x[bad] = 0
		this = Matrix::drop0(this)
		
		if ( length(gname) == ncol(x$data())) {
			print ("You have given me the correlating variables as gname - hope that was intended!")
			goi <- as.numeric(gname)
		}else {
			goi <- as.vector(t(this[gname,]))
		}
		cor.values <-  FastWilcoxTest::CorMatrix( this, goi)
		names(cor.values) = rownames(x)
		
		#calc.cor <- function(v,comp){
		#	cor(v,comp)
		#}
		#cor.values <- apply(x$data(),1,calc.cor,comp=goi)
	
		sort(cor.values)
	
} )

