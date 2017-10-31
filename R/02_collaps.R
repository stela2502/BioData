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
if ( ! isGeneric('collaps') ){ setGeneric('collaps', ## Name
		function (dataObj, by=c('median','mean','sd','sum' ), groupCol='GroupID', copy=FALSE ) { ## Argumente der generischen Funktion
			standardGeneric('collaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
}else {
	print ("Onload warn generic function 'collaps' already defined - no overloading here!")
}

setMethod('collaps', signature = c ('BioData'),
		definition = function (dataObj, by=c('median','mean','sd','sum', 'var' ), groupCol='GroupID', copy=FALSE ) {
			if ( copy ) {
				dataObj = dataObj$clone()
			}
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
			
			u <- unique(as.vector(dataObj$samples[,groupCol]))			
			m <- length(u)
			
			merged_df <- function ( orig_df ){
				mm <-  matrix ( rep(0,m * nrow(orig_df) ), ncol=m)
				colnames(mm) <- u
				rownames(mm) <- rownames(orig_df)
				new_samples <- NULL
				for ( i in u ){
					all <- which(as.vector(dataObj$samples[, groupCol]) == i )
					new_samples <- rbind ( new_samples, dataObj$samples[all[1],] )
					mm[,i] <- apply( orig_df[ , all],1,f)
				}
				list( data.frame(mm), new_samples )
			}
			
			mm <- merged_df( dataObj$dat )
			new_samples <- mm[[2]]
			mm <- mm[[1]]
			mr <- NULL
			mz <- NULL
			if (! is.null( dataObj$raw ) ){
				mr <- merged_df( dataObj$raw )
				mr <- mr[[1]]
			}
			if ( ! is.null(dataObj$zscored)) {
				mz <- merged_df( dataObj$zscored )
				mz <- mz[[1]]
			}
			
			name = paste( unlist(strsplit( paste( dataObj$name, groupCol, by, sep='_') , '\\s')) , collapse='_')
			
			try ( { reduceTo ( dataObj, what='col',
					to = setdiff( 
							as.character(dataObj$samples[,dataObj$sampleNamesCol] ), 
							as.character(new_samples[,dataObj$sampleNamesCol] )
					), 
					name= name
			) }, silent=TRUE )
			colnames(mm) <- as.vector(new_samples[, dataObj$sampleNamesCol])
			dataObj$dat <- mm
			if ( ! is.null(mr)){
				dataObj$raw = mr
			}
			if ( ! is.null(mz)){
				dataObj$zscored = mz
			}
			dataObj
})




