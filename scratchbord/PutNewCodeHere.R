SQLite_ExpressionSummary <- function (fname ) {
	
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
}

SQLite_SampleSummary <- function (fname ) {
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
}