#' @name DimReduction
#' @aliases DimReduction,BioData-method
#' @rdname DimReduction-methods
#' @docType methods
#' @description Create the initial dimension reduction dataset (PCA n =100)
#' @param x the BioData object
#' @param genes create the dim red for the genes dimension (rows) default=FALSE
#' @param n  how many eigenvectors to collapse to default=100
#' @param method which method ('auto', 'irlba', 'bpca')
#' @param force re-produced the dataset even if it already exists default=FALSE
#' @title initial dimensional reduction step based on PCA
#' @return the name of the result object in the usedObj list.
#' @export DimReduction
if ( ! isGeneric('DimReduction') ){setGeneric('DimReduction', ## Name
	function ( x, genes=FALSE, n=100, method=c('auto','irlba', 'bpca'), force=FALSE ) { 
		standardGeneric('DimReduction')
	}
) }

setMethod('DimReduction', signature = c ('BioData'),
	definition = function ( x, genes=FALSE, n=100, method=c('auto','irlba', 'bpca'), force=FALSE ) {
	
	PCA_name = 'pr'
	cmpTo = colnames(x$dat)
	if( genes) {
		PCA_name = 'prGenes'
		cmpTo = rownames(x$dat)
	}
	if ( n > length(cmpTo) ) {
		n = length(cmpTo) -1
		print ( paste("n set to",n) )
	}
	rerun = 0
	if ( is.na(match( PCA_name , names(x$usedObj))) ){
		rerun = 1
	}else{
		if (isS4(x$usedObj[[PCA_name]])) {
			## check that this object comes from the right dataset
			if ( all.equal( rownames(x$usedObj[[PCA_name]]@scores), cmpTo ) == F ) {
				rerun = 1
			}
		}else {
			if ( all.equal( rownames(x$usedObj[[PCA_name]]$x), cmpTo ) == F ) {
				rerun = 1
			}
		}
	}
	if ( force )
		rerun = 1
	
	if ( rerun == 1) {
		tmp = x$data()
		bad = which(tmp@x < 0)
		if ( length(bad) > 0 ) {
			tmp[bad] = 0
		}
		if ( ! genes ) {
			tmp = t(tmp)
		}
		if ( method == 'auto' ){
			if ( nrow(x$dat) * ncol(x$dat) > 1e6 )
				method = 'irlba'
			else 
				method= 'bpca'
		}
		if ( method == 'irlba' ) {
			message ( "irlba::prcomp_irlba is used to save memory and time (more than 1e+6 values)" )
			x$usedObj[[PCA_name]] <- irlba::prcomp_irlba ( tmp, center=T, n=n )
			rownames(x$usedObj[[PCA_name]]$x) = cmpTo
			
		}else if ( method == 'bpca') {
			message ( "pcaMethods::bpca to also use the holes (0 values)" )
			x$usedObj[[PCA_name]] <- pcaMethods::bpca( as.matrix(tmp), nPcs=n )
			rownames(x$usedObj[[PCA_name]]@scores) = cmpTo
		}else {
			stop( paste("method", method, "is not defined here" ) )
		}
		rm(tmp)
		gc()
	}
	PCA_name
} )
