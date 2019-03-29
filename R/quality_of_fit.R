#' @name quality_of_fit
#' @aliases quality_of_fit,BioData-method
#' @rdname quality_of_fit-methods
#' @docType methods
#' @description Calculates a quality of fit
#' @param obj The BioData object
#' @param what cells or genes - which grouping do you want to check?
#' @param col group name as one column in the samples table
#' @title description of function quality_of_fit
#' @export 
if ( ! isGeneric('quality_of_fit') ){ methods::setGeneric('quality_of_fit', ## Name
			function ( x, col, what='cells' ) { 
				standardGeneric('quality_of_fit')
			}
	)
}else {
	print ("Onload warn generic function 'quality_of_fit' already defined - no overloading here!")
}


setMethod('quality_of_fit', signature = c ('BioData'),
		definition = function ( x, col, what='cells' ) {			
			if ( x$logged) {
				difference <- function ( x, clusters, groups ) {
					ret = 0
					x <- exp(x)
					for ( i in groups  ) {
						a <- x[which( clusters == i)]
						a <- a[- (is.na(a))==F]
						if ( length(a) > 1 ) {  
							
							ret = ret + sum( (a- mean(a) )^2 ) 
						}
					}
					ret
				}
			}else {
				difference <- function ( x, clusters, groups.n ) {
					ret = 0 
					for ( i in groups.n  ) {
						a <- x[which( clusters == i)]
						a <- a[- (is.na(a))==F]
						if ( length(a) > 1 ) {  ret = ret + sum( (a- mean(a) )^2 ) }
					}
					ret
				}
			}
			test <- as.matrix(x$data())
			rem <- which(test ==  -20 )
			if ( length(rem) == 0) {
				rem <- which(test ==  0 )
			}
			test[ rem ] = NA
			if ( what=='cells') {
				clusters <- x$samples[,col]
				ret <- list ( 'single' = apply(test,2, difference, clusters, unique(as.character(clusters)) ) )
				ret$sum = log( round(sum(ret$single)) )
			}
			else if ( what=='genes') {
				clusters <- x$annotation[,col]
				ret <- list ( 'single' = apply(test,1, difference, as.character(clusters), max(clusters) ) )
				ret$sum = round(sum(ret$single))
			}
			else {
				stop(paste( what,': what option is not supported!') )
			}
			ret
		} 
)


