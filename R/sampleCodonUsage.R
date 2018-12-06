#' @name sampleCodonUsage
#' @aliases sampleCodonUsage,tRNAMINT-method
#' @rdname sampleCodonUsage-methods
#' @docType methods
#' @description Calculates the sum of tRNA fragments for one sample that come from a specific tRNA for one codon.
#' @param x the tRNAMINT object
#' @param sname the sample you want to depict
#' @param codons the codon list to use default=all
#' @param min_reads use all tRNA fragments with at least that many reads for the sample default=1
#' @param tRF.reliability MINT reports an exclusive and ambiguous set of tRNA fragments state one of them or use both default=NULL == both
#' @param tRF.type MINT reports different types of tRNA fragments (3'-half, 5'-half, ...) default=NULL == use all
#' @title description of function sampleCodonUsage
#' @export 
if ( ! isGeneric('sampleCodonUsage') ){ setGeneric('sampleCodonUsage', ## Name
	function ( x, sname, codons=NULL, min_reads=1, tRF.reliability=NULL, tRF.type=NULL ) { 
		standardGeneric('sampleCodonUsage')
	}
)
}else {
	print ("Onload warn generic function 'sampleCodonUsage' already defined - no overloading here!")
}

setMethod('sampleCodonUsage', signature = c ('tRNAMINT'),
	definition = function ( x, sname, codons=NULL, min_reads=1, tRF.reliability=NULL, tRF.type=NULL ) {
	
	if ( !is.null( tRF.type ) ){
		u <- grep ( tRF.type, x$annotation$tRF.type.s. )
	}else{
		u <- 1:nrow(x$annotation)
	}
	if ( !is.null( tRF.reliability ) ){
		u <- intersect ( u , grep ( tRF.reliability, x$annotation$reliability ) )
	}
	if ( is.null(codons) ) {
		codons = x$usedObj$Codons
	}
	
	u <- intersect ( u, which(x$data()[,sname] >= min_reads) )
	
	if ( length(u) == 0 ) {
		stop ("Sorry, but I could not find any read that matches your specification" )
	}
	
	apply(x$annotation[u,codons],2,sum)
} )
