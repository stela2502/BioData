#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @description create a NGSexpressionSet from a counts list object
#' @param dat the counts list you get from Rsubread::featureCounts()
#' @title description of function as_BioData
#' @export 
setGeneric('as_BioData', ## Name
	function ( dat ) { ## Argumente der generischen Funktion
		standardGeneric('as_BioData') ## der Aufruf von standardGeneric sorgt für das_BioData Dispatching
	}
)

setMethod('as_BioData', signature = c ('list'),
	definition = function ( dat ) {
	ret = NULL
	if (all.equal( names ( dat), c("counts" ,"annotation", "targets", "stat")  ) ) {
		samples <- data.frame(t(dat$stat))
		colnames(samples) <- as.character(t(samples[1,]))
		samples <- samples[-1,]
		samples$filename <- rownames(samples)
		rownames(samples) <-1:nrow(samples)
		ret <- BioData$new( 
				dat= cbind(dat$annotation, dat$counts), 
				Samples = samples, 
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
