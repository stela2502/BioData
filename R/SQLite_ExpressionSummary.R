#' @name SQLite_ExpressionSummary
#' @aliases SQLite_ExpressionSummary,BioData-method
#' @rdname SQLite_ExpressionSummary-methods
#' @docType methods
#' @description Connect to a sqlite database and retrieve a gene read count summary for all the stored genes.
#' @param fname the file to the sqlite db
#' @title description of function SQLite_ExpressionSummary
#' @export 
setGeneric('SQLite_ExpressionSummary', ## Name
	function (fname ) { ## Argumente der generischen Funktion
		standardGeneric('SQLite_ExpressionSummary') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('SQLite_ExpressionSummary', signature = c ('character'),
	definition = function (fname ) {
	
	dbh <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=fname )
	sth <- RSQLite::dbSendQuery(dbh, paste(  
					"SELECT gene_id , avg( value), count(value), gname" ,
					"from  datavalues left join genes on gene_id = genes.id",
					#			"where sample_id IN (select id from samples where sname not like '%spliced%')",  
					"GROUP by gene_id"
			)
	)
	ret <- RSQLite::dbFetch(sth)
	ret
} )
