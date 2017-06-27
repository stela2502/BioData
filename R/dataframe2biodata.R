#' @name dataframe2biodata
#' @aliases dataframe2biodata,BioData-method
#' @rdname dataframe2biodata-methods
#' @docType methods
#' @description This function creates a BioData object from a data frame only assuming, that the first column of the data frame
#' contains the row names. The sample table will be constructed from the colnames(x)[-1] vector.
#' @param x the data table
#' @param name the name of the new BioData object
#' @title description of function dataframe2biodata
#' @export 
setGeneric('dataframe2biodata', ## Name
	function (x, name="BioData") { ## Argumente der generischen Funktion
		standardGeneric('dataframe2biodata') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('dataframe2biodata', signature = c ('data.frame'),
	definition = function (x, name="BioData") {
	Samples <- data.frame( SampleName = colnames(x)[-1] )
	ret <- BioData$new( dat=x, Sample=Samples, namecol='SampleName', namerow=  colnames(x)[1], outpath=pwd(), name=name)
	ret
} )
