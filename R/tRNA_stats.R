#' @name tRNA_stats
#' @aliases tRNA_stats,tRNAMINT-method
#' @rdname tRNA_stats-methods
#' @docType methods
#' @description 
#' @param x  TEXT MISSING
#' @param acol  TEXT MISSING
#' @param scol  TEXT MISSING
#' @param norm.type  TEXT MISSING default=NULL
#' @param codon  TEXT MISSING default=NULL
#' @param fun  TEXT MISSING default=function(x) { sum(x, na.rm=TRUE) }
#' @title description of function tRNA_stats
#' @export 
setGeneric('tRNA_stats', ## Name
	function ( x, acol, scol, norm.type=NULL, codon=NULL, fun=function(x) { sum(x, na.rm=TRUE) } ) { ## Argumente der generischen Funktion
		standardGeneric('tRNA_stats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('tRNA_stats', signature = c ('tRNAMINT'),
	definition = function ( x, acol, scol, norm.type=NULL, codon=NULL, fun=function(x) { sum(x, na.rm=TRUE) } ) {
		ret <- NULL
	
	b <- x$clone()
	if ( length(unique(x$samples[,scol])) != 2){
		stop( "Sorry you can only get this statistics on exactly two sample groups")
	}
	if ( !is.null(norm.type) ) {
		reduceTo(b,what='col', colnames(x$data)[ grep( norm.type, b$samples$NormalizationMode ) ] )
	}else {
		norm.type= "All"
	}
	if ( ! is.null(codon) ) {
		reduceTo(b,what='row', rownames(b$data)[ which( b$annotation[,codon] == 1 ) ] )
	}else {
		codon = "all codons"
	}
	
	for ( name in unique(b$annotation[,acol])) {
		a <- b$clone()
		reduceTo( a, what='row', to= rownames(a$data)[which(a$annotation[,acol] == name)])
		collapse(a,what='col',group=scol, fun = function(x){ x[is.na(x)] = 0; mean(x)} )
		## now I have two columns of mean values for one state in the table
		## therefore I can add exactly one stat entry
		if ( class(a$data) == 'numeric' || nrow(a$data) < 2){
			ret <- rbind( ret, c(name, rep( NA, 7) ) )
		}else {
			print ( a$data )
			t <- wtd.t.test( x= a$data[,1], y=a$data[,2], weight= apply(a$data, 1,sum) ,samedata=TRUE )
			ret <- rbind( ret, c( name, t$coefficients,t$additional ) )
		}
	}
	colnames(ret) <- c( acol, 't.value', 'df' ,'p.value', 'Difference', paste('Mean.',a$samples[1,scol],sep=""), paste('Mean.',a$samples[2,scol],sep=""), 'Std. Err')
	if ( is.null(x$usedObj$tRNA_stats) ) {
		x$usedObj$tRNA_stats <- list ()
	}
	x$usedObj$tRNA_stats[[paste( codon, norm.type)]] <- ret
	invisible(x)
} )
