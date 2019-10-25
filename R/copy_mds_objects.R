#' @name copy_mds_objects
#' @aliases copy_mds_objects,BioData-method
#' @rdname copy_mds_objects-methods
#' @docType methods
#' @description copy MDS struxtures from one BioData object to another.
#' By default the name of the otehr BioData object is combined to the MDS name
#' @param x the BioData object the MDS tables should be copied to
#' @param from the BioData object that provides the MDS tables 
#' @param name name the MDS you want to copy (default NULL copy all)
#' @param nameExt the attached name (default from$name)
#' @title copy MDS structures and checking that the dimensions and order is OK.
#' @export 
setGeneric('copy_mds_objects', ## Name
	function ( x, from, name=NULL, nameExt=NULL) { ## Argumente der generischen Funktion
		standardGeneric('copy_mds_objects') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('copy_mds_objects', signature = c ('BioData'),
	definition = function ( x, from, name=NULL, nameExt=NULL) {
	if ( is.null(nameExt)){
		nameExt = from$name
	}
	useOnly = match( colnames(x), colnames(from) )
	if ( all(unlist(lapply(useOnly, is.na))) ){
		stop( paste( n,"->", n2, "sample names do not overlap between", x$name, "and", from$name ))
	}
	if ( length(which(is.na(useOnly))) >0){
		stop( "Warning: MDS objects do not overlap 100% (missing entries in from)")
	}		
	
	for ( listName in names(from$usedObj)[ grep( '^MDS', names(from$usedObj)) ] ) {
		for (n in names(from$usedObj[[listName]])){
			if ( !is.null( name )) {
				if ( n != name ) {
					next;
				}				
			}
			if ( is.null(x$usedObj[[listName]])){ x$usedObj[[listName]] = list()}
			new_name = n
			if ( nameExt != "" ) {
				new_name = paste( nameExt, n, sep="_" )
			}
			x$usedObj[[listName]][[new_name]] = from$usedObj[[listName]][[n]][useOnly,]
		}
	}
	 
	invisible(x)
} )
