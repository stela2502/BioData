#' To identify variable genes is almost the first thing that every single cell analysis strategy does.
#' The existing VarGene function identify many highly expressed genes like 
#' ribosomal genes or mitochondrial genes if they have not been removed before.
#' 
#' This normally produces problems in the downstream analysis as they normally correlate with the cell cycle which
#' most of the times is not the main interesting cellular process.
#' 
#' Here we try to identify genes that are correlated less to the UMI count than would be expected.
#' The function estimates the normal correlation values for each gene count and selects for genes that
#' are correlated to the UMI count less than other genes with the same fraction of cells expressing them.
#' 
#' This has proven to especially pick genes at the lower expression range like 
#' chromatin remodeling proteins or transcription factors.
#' 
#' No normalization is required for this function to work.
#' It is always applied to the not normalized data.
#' 
#' @name getGenesExpressedHigherThanExpected
#' @aliases getGenesExpressedHigherThanExpected,BioData-method
#' @rdname getGenesExpressedHigherThanExpected-methods
#' @docType methods
#' @description My personal VarGenes function.
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
				
		OrderN = paste('GExHiEx_OoI',name,sep="_")
		DataN = paste('GExHiEx_RES',name,sep="_")
		
		if ( sum(is.na(match( c(OrderN, DataN), colnames(x$annotation) )))==0 ){	
		#if ( ! is.null(x$annotation[, OrderN]) & ! is.null(x$annotation[, DataN])){
			print( paste( "The columns", DataN, "and", OrderN,"are existsing - no re-calculation") )
			return( sort(rownames(x)[which(x$annotation[,OrderN] <= n )]))
		}
		
		obj <- reduceTo( x, copy=T, what='col', to= colnames(x)[which(x$samples[,group] == name)], name=name)
		if ( is.null(x$raw)){
			mat = obj$dat
		}else {
			mat = obj$raw
		}
		cor = FastWilcoxTest::CorMatrix( mat, obj$samples[, ref] )
		names(cor) = rownames(obj)
		B = FastWilcoxTest::ColNotZero( Matrix::t(mat))
		g = lm( cor ~ B )
		# get the data
		order_of_interest = names(sort(g$residuals))
		OrderN = paste('GExHiEx_OoI',name,sep="_")
		DataN = paste('GExHiEx_RES',name,sep="_")
		# init order vector
		x$annotation[, OrderN] = rep(nrow(x$dat), nrow(x$dat))
		# init data vector
		x$annotation[, DataN] = rep(0, nrow(x$dat))
		
		m = match( order_of_interest, rownames(x) )
		x$annotation[m, OrderN] = 1:length(m)
		x$annotation[m, DataN ] = order_of_interest
		sort(order_of_interest[1:n+1])
		
	})
	names(genes) = levels( x$samples[,group])
	genes

} )
