#' @name simpleAnova
#' @aliases simpleAnova,BioData-method
#' @rdname simpleAnova-methods
#' @docType methods
#' @description  This function calculates an annova to identify significant changes in the BioData
#' @description  has a higher sensitivity for multi group analyses to identify group specific changes
#' @description  or general trends in the dataset. This function adds the results into the stats slot
#' @description  of the BioData object.
#' @param x the BioData object
#' @param groupCol the samples table column that contains the grouping information
#' @param padjMethod the p value correction method as described in  \code{\link[stats]{p.adjust}}
#' @title description of function simpleAnova
#' @export 
setGeneric('simpleAnova', ## Name
		function ( x, groupCol='GroupName', padjMethod='BH' ) { 
			standardGeneric('simpleAnova')
		}
)

setMethod('simpleAnova', signature = c ( 'BioData') ,
		definition = function ( x, groupCol='GroupName', padjMethod='BH' ) {
			x <- normalize(x)
			significants <- apply ( x$data() ,1, function(a) { anova( lm (a ~ x$samples[,groupCol ]))$"Pr(>F)"[1] } )
			adj.p <- p.adjust( significants, method = padjMethod)
			res <- data.frame( genes= rownames(x$dat), pvalue= significants,  adj.p )
			colnames(res)[3] <- paste('padj',padjMethod)
			if ( length (x$stats) == 0 ){
				x$stats <- list ( 'simpleAnova' = res )
			}
			else {
				x$stats[[length(x$stats)+1]] <- res
				names(x$stats)[length(x$stats)] = 'simpleAnova'
			}
			x
		}
)