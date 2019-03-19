#' @name copy_mds_objects
#' @aliases copy_mds_objects,BioData-method
#' @rdname copy_mds_objects-methods
#' @docType methods
#' @description copy MDS struxtures from one BioData object to another.
#' By default the name of the otehr BioData object is combined to the MDS name
#' @param x the BioData object the MDS tables should be copied to
#' @param from the BioData object that provides the MDS tables 
#' @param nameExt the attached name (default from$name)
#' @title copy MDS structures and checking that the dimensions and order is OK.
#' @export 
setGeneric('copy_mds_objects', ## Name
	function ( x, from, nameExt=NULL) { ## Argumente der generischen Funktion
		standardGeneric('copy_mds_objects') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('copy_mds_objects', signature = c ('BioData'),
	definition = function ( x, from, nameExt=NULL) {
	if ( is.null(nameExt)){
		nameExt = from$name
	}
	lapply( names(from$usedObj)[grep ( 'MDS', names(from$usedObj))],
			function(n) {
				new_name = paste( nameExt, n, sep="_" )
				if ( is.null(x$usedObj[[n]] )){
					x$usedObj[[n]] = list()
				}
				x$usedObj[[n]] = from$usedObj[[n]]
				## fix the ordering
				lapply( names( from$usedObj[[n]] ), function(n2) {
							m = match( colnames(x$dat), rownames (from$usedObj[[n]][[n2]]))
							if ( all.equal(  colnames(x$dat), rownames (from$usedObj[[n]][[n2]])[m] ) == FALSE ) {
								stop( paste( n,"->", n2, "sample names do not overlap between", x$name, "and", from$name ))
							}
							x$usedObj[[n]][[paste(nameExt, n2, sep="_" )]] = from$usedObj[[n]][[n2]][m,]
						})      
				NULL
			})   
	invisible(x)
} )
