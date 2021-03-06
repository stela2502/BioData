#' @name collaps
#' @aliases collaps,BioData-method
#' @rdname collaps-methods
#' @docType methods
#' @description  This function will collpase the data in the BioData to only contain one value
#' per sample group.
#' @param dataObj the BioData
#' @param by collapsing method c('median','mean','sd','sum', or own function )
#' @param groupCol the sample names you want to group on
#' @param copy create a copy of the data (default=F)
#' @title description of function collaps
#' @export 
if ( ! isGeneric('collaps') ){ methods::setGeneric('collaps', ## Name
		function (dataObj, by=c('median','mean','sd','sum' ), groupCol='GroupID', copy=FALSE ) { 
			standardGeneric('collaps')
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
				by = "user_function"
			}else {
			switch( by,
					median = f<- function (x ) { median(x) },
					mean = f <- function(x) { mean(x) },
					sd = f <- function(x) { stats::sd(x) },
					sum = f <-function(x) { sum(x[which(x > 0)])},
					var = f <- function(x) { stats::var(x) }
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
					if ( length(all) == 1 ) {
						mm[,i] <- orig_df[ , all]
					}else if (  length(all) > 1) {
						mm[,i] <- apply( orig_df[ , all],1,f)
					}else {
						stop( paste( "The entry", i, "has no column in the data - please fix that!"))
					}
					
				}
				list( df = data.frame(mm), new_sample = new_samples )
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
					to = as.character(new_samples[,dataObj$sampleNamesCol] ), 
					name= name
			) }, silent=TRUE )
			dataObj$samples = new_samples	
			colnames(mm) <- as.vector(new_samples[, dataObj$sampleNamesCol])
			dataObj$dat <- Matrix::Matrix(as.matrix(mm), sparse=T)
			if ( ! is.null(mr)){
				colnames(mr) <- as.vector(new_samples[, dataObj$sampleNamesCol])
				dataObj$raw = mr
			}
			if ( ! is.null(mz)){
				colnames(mz) <- as.vector(new_samples[, dataObj$sampleNamesCol])
				dataObj$zscored = mz
			}
			changeNames(dataObj, what='col',colname= groupCol )
			invisible(dataObj)
})




