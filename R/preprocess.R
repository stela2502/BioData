#' @name preprocess
#' @aliases preprocess,BioData-method
#' @rdname preprocess-methods
#' @docType methods
#' @description create the count dataset and normalize it following DEseq standard
#' @description  function used by DEseq Do not use
#' @param x the BioData
#' @param condition the column in the samples table to use as grouing variable
#' 		  The condition will be used as formular paste('~', condition) 
#'        and can therfore contain more complex function like 'batch + condition' 
#' @title description of function preprocess
#' @export
if ( ! isGeneric('preprocess') ){ methods::setGeneric('preprocess', ## Name
		function (x, condition, add=NULL ) { 
			standardGeneric('preprocess')
		}
)
}else {
	print ("Onload warn generic function 'preprocess' already defined - no overloading here!")
}

setMethod('preprocess', signature = c ('BioData'),
		definition = function (x, condition) {
			if ( is.null( x$usedObj[['cds']]) ) {
				x$usedObj[['cds']] <- list()
			}
			
			if ( is.na( match( condition, names(x$usedObj[['cds']]) ) )) {
				#condition <- as.factor(x$samples[,condition])
				dat <- list()
				if ( nrow(x$raw) > 0 ) {
					t <- as.matrix(x$raw)
				}else {
					t <- as.matrix(x$dat)
				}
				#colnames( t ) <- x$samples[,condition]
				colnames( t ) <- make.names(colnames( t ))
				s <- x$samples
				rownames(s) <- make.names(rownames(s))
				#browser()
				#storeage.mode(t) <- 'numeric'
				
				dat$cds <- DESeq2::DESeqDataSetFromMatrix(
					t,
					x$samples,
					formula( paste('~',condition ) )
				)
				#dat$cds <- DESeq2::estimateSizeFactors(dat$cds)
				#DESeq2::sizeFactors(x$cds)
				#dat$cds <- DESeq2::estimateDispersions(dat$cds)
				#dat$vsdFull = DESeq2::varianceStabilizingTransformation( dat$cds )
				x$usedObj[['cds']][[length(x$usedObj[['cds']])+1]] = dat$cds
				names(x$usedObj[['cds']])[length(x$usedObj[['cds']])] = condition
			}
			x
		}
)
