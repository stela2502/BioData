#' @name SQLite_SampleSummary
#' @aliases SQLite_SampleSummary,BioData-method
#' @rdname SQLite_SampleSummary-methods
#' @docType methods
#' @description Connect to a sqlite database and retrieve a sample read count summary for all the stored genes.
#' @param fname the file to the sqlite db
#' @param cells details about the cell storage table default list( 'table' = 'samples', 'rev' = 'sample_id', 'name' = 'sname')
#' @title description of function SQLite_SampleSummary
#' @export 
setGeneric('SQLite_SampleSummary', ## Name
	function (fname , cells= list( 'table' = 'samples', 'rev' = 'sample_id', 'name' = 'sname') ) { 
		standardGeneric('SQLite_SampleSummary')
	}
)

setMethod('SQLite_SampleSummary', signature = c ('character'),
	definition = function (fname, cells=  list( 'table' = 'samples', 'rev' = 'sample_id', 'name' = 'sname') ) {
	dbh <- RSQLite::dbConnect(RSQLite::SQLite(),dbname=fname )
	sth <- RSQLite::dbSendQuery(dbh, paste(  
					paste("SELECT",cells$rev,", sum(value) as reads, count(value) as count,", cells$name ),
					paste("from  datavalues left join ", cells$table,"on", cells$rev, " = ",paste( sep=".", cells$table, 'id') ),
					#"where sample_id IN (select id from samples where sname not like '%spliced%')",  
					paste("GROUP by ", cells$rev)
			)
	)
	ret <- RSQLite::dbFetch(sth)
	ret
} )
