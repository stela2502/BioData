#' @name as.BioData
#' @alias.BioDataes as.BioData,BioData-method
#' @rdname as.BioData-methods
#' @docType methods
#' @description create a NGSexpressionSet from a counts list object
#' @param dat the counts list you get from Rsubread::featureCounts()
#' @title description of function as.BioData
#' @export 
setGeneric('as.BioData', ## Name
	function ( dat ) { ## Argumente der generischen Funktion
		standardGeneric('as.BioData') ## der Aufruf von standardGeneric sorgt für das.BioData Dispatching
	}
)

setMethod('as.BioData', signature = c ('list'),
	definition = function ( dat ) {
	ret = NULL
	if (all.equal( names ( dat), c("counts" ,"annotation", "targets", "stat")  ) ) {
		samples <- data.frame(t(dat$stat))
		colnames(samples) <- as.BioData.character(t(samples[1,]))
		samples$filename <- rownames(samples)
		rownames(samples) <-1:nrow(samples)
		ret <- BioData$new( 
				dat= cbind(dat$annotation, dat$counts), 
				samples = samples, 
				namecol= 'filename', 
				namerow= 'GeneID',
				outpath= ''
		)
	}
	else {
		print ("The list needs to contain the entries counts ,annotation, targets and stat" )
	}
	ret
} )
