#' @name IdentifyMarkerGenes
#' @aliases IdentifyMarkerGenes,BioData-method
#' @rdname IdentifyMarkerGenes-methods
#' @docType methods
#' @description This method uses whichever stats function was selected for this class using the cresteStats function.
#' But it compares one group versus all other groups to find marker genes for this group only.
#' This function will add length(group) new stats tables.
#' @param x the BioData object
#' @param gname the samples column name to group on.
#' @title description of function IdentifyMarkerGenes
#' @export 
setGeneric('IdentifyMarkerGenes', ## Name
	function ( x, gname ) { ## Argumente der generischen Funktion
		standardGeneric('IdenifyMarkerGenes') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('IdentifyMarkerGenes', signature = c ('BioData'),
	definition = function ( x, gname ) {
		tmp <- x$clone()
		for ( n in unique(tmp$samples[,gname])) {
			tmp$stats <- NULL
			gc(FALSE)
			new_g <-  paste( 'IdentifyMarkerGenes',  gname, n )
			print (paste( "Processing:", new_g ))
			g <- rep('rest', ncol(x$dat) )
			g[which(tmp$samples[,gname] == n )] = n
			tmp$samples[,new_g] <- factor( g, levels= c( n, 'rest') ) 
			createStats( tmp, new_g)
			x$stats[[new_g]] = tmp$stats[[1]]
		}
		rm(tmp)
		gc(FALSE)
		invisible(x)
} )
