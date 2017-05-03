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
setGeneric('collapse2codons', ## Name
	function ( x, sum.fun=function(x) { sum (x, na.rm=T)} , name="collapsed") { ## Argumente der generischen Funktion
		standardGeneric('collapse2codons') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('collapse2codons', signature = c ('tRNAMINT'),
	definition = function ( x, sum.fun=function(x) { sum (x, na.rm=T)} , name="collapsed") {
	
	dat = NULL
	for ( codon in x$usedObj$Codons ) {
		a <- x$clone()
		reduceTo(a,'row', rownames(a$data)[which(a$annotation[,codon] == 1)] )
		if ( class(a$data) == 'numeric') {
			dat <- rbind( dat, a$data )
			a$data <- t(data.frame(a$data))
		}else {
			dat <- rbind( dat, apply( a$data,2,sum.fun))
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
