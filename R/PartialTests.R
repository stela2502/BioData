#' PartialTests should take two column names and apply stat tests on sub-samples using group1 
#' and test for differentials in group2.
#' The stat results are saved inside the original R6 object
#' All differential genes are returned as list.
#' @name PartialTests
#' @aliases PartialTests,BioData-method
#' @rdname PartialTests-methods
#' @docType methods
#' @description Run tests firs subselecting the BioData object on the outer grouping and calculating tests on the inner grouping.
#' @param obj The BioData object
#' @param groupA the outer grouping column name
#' @param groupB the inner grouping column name
#' @param pcut p cutoff value for gene selection default=1e-5
#' @param logfc.threshold Cpp test option default= .1
#' @param minPct Cpp test option default= .1
#' @title description of function PartialTests
#' @export 
setGeneric('PartialTests', ## Name
	function (obj, groupA, groupB, pcut=1e-5, logfc.threshold= .1, minPct= .1 ) { ## Argumente der generischen Funktion
		standardGeneric('PartialTests') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('PartialTests', signature = c ('BioData'),
	definition = function (obj, groupA, groupB, pcut=1e-5, logfc.threshold= .1, minPct= .1 ) {
        thisN = paste('PartialTests', obj$name, groupA,groupB, logfc.threshold, minPct , sep="_")
#browser()
        if ( !is.null( obj$usedObj[[thisN]] ) ){
                ret = lapply( 
                        obj$usedObj[[thisN]], 
                        function(x) { unique( x[which(x[,'p_val_adj'] < pcut),'gene']) } 
                )
                return ( ret )

        }
        stats <- lapply( levels(obj$samples[,groupA]), 
          function( Ga ) {
            Test <- reduceTo( obj, copy=T, what='col', 
                to=colnames(obj)[which(obj$samples[, groupA] == Ga)], 
                name=paste( sep="_", obj$name, groupA, Ga) )
            if ( length( grep( groupB, names(Test$stats))) > 0 ){
                Test$stats=list()
            }
            Cpp_FindAllMarkers( Test, groupB, logfc.threshold= logfc.threshold, minPct= minPct)
            Test$stats[[  grep( groupB, names(Test$stats)) ]]
        } )
        names(stats) = paste( groupA, levels(obj$samples[,groupA]), sep=".")
        obj$usedObj[[thisN]] = stats
        print (paste('All stats have been storerd in obj$usedObj[[',thisN,']].'))
        ret = lapply( obj$usedObj[[thisN]], function(x) { unique( x[which(x[,'p_val_adj'] < pcut),'gene']) } )
        return ( ret )
} )
