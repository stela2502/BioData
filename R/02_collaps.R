#' @name collaps
#' @aliases collaps,BioData-method
#' @rdname collaps-methods
#' @docType methods
#' @description  This function will collpase the data in the BioData to only contain one value
#' @description  per sample group.
#' @param dataObj the BioData
#' @param by collapsing method c('median','mean','sd','sum', or own function )
#' @param groupCol the sample names you want to group on
#' @param copy create a copy of the data (default=F)
#' @title description of function collaps
#' @export 
setGeneric('collaps', ## Name
		function (dataObj, by=c('median','mean','sd','sum' ), groupCol='GroupID', copy=FALSE ) { ## Argumente der generischen Funktion
			standardGeneric('collaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('collaps', signature = c ('BioData'),
		definition = function (dataObj, by=c('median','mean','sd','sum', 'var' ), groupCol='GroupID', copy=FALSE ) {
			if ( copy ) {
				dataObj = dataObj$clone()
			}
			u <- unique(as.vector(dataObj$samples[,groupCol]))
			m <- length(u)
			mm <-  matrix ( rep(0,m * nrow(dataObj$data)), ncol=m)
			colnames(mm) <- u
			rownames(mm) <- rownames(dataObj$data)
			f <- NULL
			if ( is.function(by)){
				f <- by
			}else {
			switch( by,
					median = f<- function (x ) { median(x) },
					mean = f <- function(x) { mean(x) },
					sd = f <- function(x) { sd(x) },
					sum = f <-function(x) { sum(x)},
					var = f <- function(x) { var(x) }
			);
			}
			if ( is.null(f) ) {
				stop("Please set what to one of 'median','mean','sd','sum'" )
			}
			new_samples <- NULL
			for ( i in u ){
				all <- which(as.vector(dataObj$samples[, groupCol]) == i )
				new_samples <- rbind ( new_samples, dataObj$samples[all[1],] )
				mm[,i] <- apply( dataObj$data[ , all],1,f)
			}
			name = paste( unlist(strsplit( paste( dataObj$name, groupCol, by, sep='_') , '\\s')) , collapse='_')
			try ( { reduceTo ( dataObj, what='col',
					to = setdiff( 
							dataObj$samples[,dataObj$sampleNamesCol] , 
							new_samples[,dataObj$sampleNamesCol] 
					), 
					name= name
			) }, silent=TRUE )
			colnames(mm) <- as.vector(new_samples[, dataObj$sampleNamesCol])
			dataObj$data <- data.frame(mm)
			dataObj
})



