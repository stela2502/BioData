#' @name stat_table
#' @aliases stat_table,tRNAMINT-method
#' @rdname stat_table-methods
#' @docType methods
#' @description converts a list of tables one table adding the list name as first column
#' @param l the list to convert
#' @param padjMethod the method name for the p.adjust call (default = 'BH')
#' @title description of function stat_table
#' @export 
setGeneric('stat_table', ## Name
	function ( l, padjMethod='BH' ) { ## Argumente der generischen Funktion
		standardGeneric('stat_table') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('stat_table', signature = c ('list'),
	definition = function ( l, padjMethod='BH' ) {
	ret <- NULL
	for ( name in names(l)){
		if ( class(l[[name]]) == 'matrix'){
			ret <- rbind( ret, cbind( rep(name, nrow(l[[name]])), l[[name]]))
		}     
	}
	colnames(ret)[1] <- 'Codon'
	ret <- data.frame(ret)
	if ( ! is.null( match('p.value', ret) ) ){
		ret <- cbind( ret, 'p.adj' =p.adjust( as.numeric(as.character(ret$p.value)), method = padjMethod) )
	}
	ret   
} )
