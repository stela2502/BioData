#' @name collapse2codons
#' @aliases collapse2codons,tRNAMINT-method
#' @rdname collapse2codons-methods
#' @docType methods
#' @description This function collapses the data to the different codons and reports the sum.fun of all codons.
#' Single fragments can/will be counted multiple times.
#' @param x the tRNAMINT object
#' @param sum.fun how to sum up the values for one sample default=function(x) { sum (x, na.rm=T)} 
#' @param name the name of the collapsed object default="collapsed"
#' @title description of function collapse2codons
#' @export 
if ( ! isGeneric('collapse2codons') ){ setGeneric('collapse2codons', ## Name
	function ( x, sum.fun=function(x) { sum (x, na.rm=T)} , name="collapsed") { ## Argumente der generischen Funktion
		standardGeneric('collapse2codons') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'collapse2codons' already defined - no overloading here!")
}

setMethod('collapse2codons', signature = c ('tRNAMINT'),
	definition = function ( x, sum.fun=function(x) { sum (x, na.rm=T)} , name="collapsed") {
	print ( "please doublecheck this function - not changed in the data to dat conversion!!!!!")
	dat = NULL
	for ( codon in x$usedObj$Codons ) {
		a <- x$clone()
		reduceTo(a,'row', rownames(a$dat)[which(a$annotation[,codon] == 1)] )
		if ( class(a$dat) == 'numeric') {
			dat <- rbind( dat, a$dat )
			a$dat <- t(data.frame(a$dat))
		}else {
			dat <- rbind( dat, apply( a$dat,2,sum.fun))
		}
	}
	dat <- data.frame(dat)
	
	ret <- tRNAMINT$new( 
			dat = cbind(Codon = x$usedObj$Codons, dat ), 
			Samples = x$samples, 
			name=name, 
			namecol=x$sampleNamesCol, 
			namerow= 'Codon', 
			outpath=x$outpath 
	)
	invisible(ret)
} )
