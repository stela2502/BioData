#' After excessive clustering and really good grouping of the cells this function can be apllied to single cells dataset
#' to identify genes that are expressed to a higher level than what the correlation to the UMI count can explain.
#' 
#' This strategy is hily questionable, but results in quite interesting gene lists. No normalization is required for this function to work.
#' It is always applied to the not normalized data.
#' 
#' @name getGenesExpressedHigherThanExpected
#' @aliases getGenesExpressedHigherThanExpected,BioData-method
#' @rdname getGenesExpressedHigherThanExpected-methods
#' @docType methods
#' @description 
#' @param x a SingleCells object
#' @param group a samples grouping used to subset the original object
#' @param ref the correlating samples group (default nUMI)
#' @param n get the top interesting genes default= 300 (in each group!)
#' @title description of function getGenesExpressedHigherThanExpected
#' @export 
if ( ! isGeneric('getGenesExpressedHigherThanExpected') ){setGeneric('getGenesExpressedHigherThanExpected', ## Name
	function ( x, group, ref='nUMI', n= 300) { 
		standardGeneric('getGenesExpressedHigherThanExpected')
	}
) }

setMethod('getGenesExpressedHigherThanExpected', signature = c ('SingleCells'),
	definition = function ( x, group, ref='nUMI', n= 300 ) {
		
	genes = lapply ( levels( x$samples[,group]) , function(name) {
		obj <- reduceTo( x, copy=T, what='col', to= colnames(x)[which(x$samples[,group] == name)], name=name)
		if ( is.null(x$raw)){
			mat = x$dat
		}else {
			mat = x$raw
		}
		cor = FastWilcoxTest::CorMatrix( mat, obj$samples[, ref] )
		names(cor) = rownames(obj)
		B = FastWilcoxTest::ColNotZero( Matrix::t(mat))
		g = lm( cor ~ B )
		SD = sd( g$residuals )
		#browser()
		order_of_interest = names(sort(g$residuals))
		x$annotation$orderOfInterest_GenesExpressedHigherThanExpected = nrow(x$dat)
		m = match( order_of_interest, rownames(x) )
		x$annotation$orderOfInterest_GenesExpressedHigherThanExpected [m] = 1:length(m)
		sort(order_of_interest[1:n])
		#sort(names(which(g$residuals < - cut * SD )))
	})
	names(genes) = levels( x$samples[,group])
	genes

} )
