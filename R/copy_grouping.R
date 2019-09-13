#' @name copy_grouping
#' @aliases copy_grouping,BioData-method
#' @rdname copy_grouping-methods
#' @docType methods
#' @description Copy a grouping from a (partial) BioData object
#' @param x the traget object
#' @param from the source object
#' @param gname the (old) group name in the source object
#' @param newgname the (new) group name for the taget object default=NULL; use gname
#' @title Copy a grouping (samples column) from one to another BioData object
#' @export 
setGeneric('copy_grouping', ## Name
	function ( x, from, gname, newgname=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('copy_grouping') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('copy_grouping', signature = c ('BioData'),
	definition = function ( x, from, gname, newgname=NULL ) {
	# add to object 'x' from object 'from' the new(?) group 'gname' and use the 'newgname' grouping name.
	if ( is.na(match(gname, colnames(from$samples)))){
		stop(paste("The group", gname, "is not defined in the from object",from$name))
	}
	if ( is.null(newgname) ) {
		newgname = gname
	}
	if ( is.na(match( newgname, colnames(x$samples) ))) {
		x$samples[,newgname] = NA
	}
	m = match( colnames(from), colnames(x) ) 
	if ( length(which(is.na(m))) > 0 ){
		print( paste("There are cells in object", from$name, "that are missing in the target object", x$name ))
		m = m[-which(is.na(m))]
	}
	if ( ! all( is.na(x$samples[m,newgname]))) {
		stop( paste("the action would delete data in the target object",x$name, " - STOP"))
	}else {
		x$samples[m,newgname] = 
				factor( paste(sep="_", from$name, from$samples[,gname] ), 
					levels = paste(sep="_", from$name, levels(from$samples[,gname]) )) 
	}
	colors_4( x, newgname, force=T)
	invisible(x)
} )
