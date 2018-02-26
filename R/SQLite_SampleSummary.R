#' @name SQLite_SampleSummary
#' @aliases SQLite_SampleSummary,BioData-method
#' @rdname SQLite_SampleSummary-methods
#' @docType methods
#' @description Connect to a sqlite database and retrieve a sample read count summary for all the stored genes.
#' @param fname the file to the sqlite db
#' @title description of function SQLite_SampleSummary
#' @export 
setGeneric('SQLite_SampleSummary', ## Name
	function (fname ) { ## Argumente der generischen Funktion
		standardGeneric('SQLite_SampleSummary') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('SQLite_SampleSummary', signature = c ('character'),
	definition = function (fname ) {
	dbh <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=fname )
	sth <- RSQLite::dbSendQuery(dbh, paste(  
					"SELECT sample_id , sum(value) as reads, count(value) as count, sname" ,
					"from  datavalues left join samples on sample_id = samples.id",
					#"where sample_id IN (select id from samples where sname not like '%spliced%')",  
					"GROUP by sample_id"
			)
	)
	ret <- RSQLite::dbFetch(sth)
	ret
} )
