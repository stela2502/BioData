#' @name tRNA_stats
#' @aliases tRNA_stats,tRNAMINT-method
#' @rdname tRNA_stats-methods
#' @docType methods
#' @description This function calculates the tRNA statistics as described in XYZ
#' @param x the tRNAMINT object
#' @param acol the statistics will be calculated for every entry in this column
#' @param scol column containing the sample grouping (two groups are required)
#' @param norm.type a tRNAMINT object normally contains 3 or more different norm datasets per sample.
#' This value is used as a patter to select one normalization mode.
#' @param codon the tRNAMINT object can extract codon information from the fragments using the extractCodonInformation() function.
#' If you give one codon column name here (x$usedObj$Codons) the stats will be calculated for this codon only
#' @param fun multiple fragments are summed up for this analysis; this is the summing function: default=function(x) { sum(x, na.rm=TRUE) }
#' @title description of function tRNA_stats
#' @export 
setGeneric('tRNA_stats', ## Name
	function ( x, acol, scol, norm.type=NULL, codon=NULL, fun=function(x) { x[is.na(x)] = 0; mean(x) } ) { ## Argumente der generischen Funktion
		standardGeneric('tRNA_stats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('tRNA_stats', signature = c ('tRNAMINT'),
	definition = function ( x, acol, scol, norm.type=NULL, codon=NULL, fun=function(x) { x[is.na(x)] = 0; mean(x) } ) {
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
	return.val <- list()
	for ( name in unique(b$annotation[,acol])) {
		a <- b$clone()
		reduceTo( a, what='row', to= rownames(a$data)[which(a$annotation[,acol] == name)])
		if ( class(a$data) !=  "matrix" & class(a$data) !=  "data.frame" ) {
			next ## useless as we have less than 2 values for this
		}
		collapse(a,what='col',group=scol, fun = fun )
		a$name = paste(codon, name)
		return.val[[paste(codon, name)]] <- a
		## now I have two columns of mean values for one state in the table
		## therefore I can add exactly one stat entry
		if ( class(a$data) == 'numeric' || nrow(a$data) < 2){
			ret <- rbind( ret, c(name, rep( NA, 8) ) )
		}else {
			if ( nrow(a$data) > 2 ) {
				#t <- cor.test( as.vector(a$data[,1]),as.vector(a$data[,2])) # this is a naive test! useless!
				#t <- wtd.t.test( x= a$data[,1], y=a$data[,2], weight= apply(a$data, 1,sum) ,samedata=TRUE )
				t <- t.test(x= a$data[,1], y=a$data[,2], paired=T)
				t$foldChange <- mean(apply(a$data,1,function(x) { x[1] / x[2]} ) )
				#ret <- rbind( ret, c( name, nrow(a$data), as.numeric(t(unlist(t))[1:4]) ) )
				#ret <- rbind( ret, c( name, t$coefficients,t$additional ) )
				ret <- rbind( ret, c( name, nrow(a$data), t$statistic,t$parameter, t$p.value, t$alternative, t$method, t$estimate, t$foldChange) )
				
			}else {
				ret <- rbind( ret, c( name, nrow(a$data), rep(NA,7) ) )
			}
		}
	}
	if ( class(ret) != "matrix") {
		print (paste("Not enough data for the fragment types ", paste( unique(b$annotation[,acol]) , collapse=", "), "and codon", codon ))
	}else {
		try({
					#colnames(ret) <- c( acol, 'frament [n]', 'statistic.t', 'df' ,'p.value', 'estimate.cor' )
					#colnames(ret) <- c( acol, 't.value', 'df' ,'p.value', 'Difference', 
					#	paste('Mean.',a$samples[1,scol],sep=""), 
					#	paste('Mean.',a$samples[2,scol],sep=""), 'Std. Err'
					#			)
					colnames(ret) <- c(acol, 'frament [n]', "t.statistic","t.parameter", "p.value", "t.alternative", "t.method", "t.estimate", 'manual.foldChange' )
					if ( is.null(x$usedObj$tRNA_stats) ) {
						x$usedObj$tRNA_stats <- list ()
					}
					x$usedObj$tRNA_stats[[paste( codon, norm.type)]] <- ret
				}
		)
	}
	invisible(return.val)
	#invisible(x)
} )
