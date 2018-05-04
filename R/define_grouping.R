#' @name define_grouping
#' @aliases define_grouping,BioData-method
#' @rdname define_grouping-methods
#' @docType methods
#' @description This will create a new sample or annotation column based on an old grouping and a list.
#' @param x The BioData object
#' @param basedOn the old grouping name
#' @param colname the new grouping name
#' @param newNames the list "with new_group_name" = c( "old Group Names" ) structure
#' @param what create a sample ('col') or annotation ('row') grouping
#' @title description of function define_grouping
#' @example 
#' # data is a BioData object with a sample group 'oldG' and the levels c( 1, 2, 3, 4)
#' define_grouping( data, 'oldG', 'newGrouping', list( 'A' = c( 1,4), B= c(2,3)), 'col')
#' # will create a 'newGrouping' column in the $samples table with the entries A and B
#' @export 
setGeneric('define_grouping', ## Name
	function ( x, basedOn, colname, newNames, what='col' ) { ## Argumente der generischen Funktion
		standardGeneric('define_grouping') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('define_grouping', signature = c ('BioData'),
	definition = function ( x, basedOn, colname, newNames, what='col' ) {
	if ( ! class(newNames) == 'list' ) {
		stop( "I need a list as newNames object" )
	}else {
		t <- table( unlist(newNames) )
		m <- which( t != 1 )
		if ( length(m) != 0 ){
			stop(
					paste( "Sorry, your new names have problematic entres:", 
							paste( 
									paste( names(t)[m],t[m], sep="="), collapse=" ,"
							) 
					)
			)
		}
	}
	
	selectFrom = NULL
	if ( what == 'col' ){
		selectFrom = x$samples
	}else if ( what == 'row') {
		selectFrom = x$annotation
	}else {
		stop( paste( "Sorry, but the what option '",what,"' is not supported", sep="" ) )
	}
	
	if ( is.null(selectFrom[,basedOn]) ) {
		stop( paste( "sample/annotation column",basedOn, "is not defined in the object"))
	}
	old_grouping <- levels(selectFrom[,basedOn])
	new <- rep ( 'noName', nrow( selectFrom ) )
	for ( n in names(newNames) ) {
		m <- match( selectFrom[,basedOn], newNames[[n]] )
		if ( length(which(!is.na(m))) > 0 ){
			m <- which(!is.na(m))
			new[m] = n
		}
	}
	if ( is.na(table(new)['noName']) ){
		if ( what == 'col' ){
			x$samples[,colname] = factor(new, levels=names(newNames))
		}else if ( what == 'row') {
			x$annotation[,colname] = new
		}
	}else {
		stop( paste( "you miss the old group names", paste( unique(t(x$samples[ grep('noName', new),basedOn] ) ), collapse=", " ), "in your newNames list"  ) )
	}
	colors_4( x, colname )
	invisible(x)
	
}  )
