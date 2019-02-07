#' @name bestGrouping
#' @aliases 'bestGrouping,BioData-method
#' @docType methods
#' @description The function is using a randomForest classifier with 2000 trees to classify the given data using the given grooping
#' @description All groups that fail to be prediceted using the random forest are deemed ungrouped.
#' @description All groups where less than 50 percent of the total samples geting classified as being from that group fail.
#' @param x the single cells ngs object
#' @param group a vector of sample columns that should be checked (the most complex is used only)
#' @param bestColname the column name to store the best grouping in
#' @param cutoff the cutoff percentage where all groups showing less than this percentacge of remapped samples are dropped
#' @param ntree for the rf process (default 2000)
#' @param ... additional variables for the randomForest call
#' @title description of function randomForest
#' @return a distRF object to be analyzed by pamNew
#' @examples
#' ## lets data be a source BioData object and selection be a smaller one with only selected samples 
#' ## for which a grouping has been produced named 'new_perfect_grouping'
#' predictor <- bestGrouping( selection, 'new_perfect_grouping' )
#' data$samples$new_perfect_grouping <- predict(predictor, as.matrix(t(data$dat)))
#' @export 
if ( ! isGeneric('bestGrouping') ){ setGeneric('bestGrouping',
		function ( x, group , bestColname='QualifiedGrouping', cutoff=0.5, ntree=2000, ...){
			standardGeneric('bestGrouping')
		}
)
}else {
	print ("Onload warn generic function 'bestGrouping' already defined - no overloading here!")
}
setMethod('bestGrouping', signature = c ('BioData'),
		definition = function (x, group, bestColname='QualifiedGrouping' , cutoff=0.5, ... ) {
			uObj <- paste( 'Predictive RFobj', group )
			rf <- NULL
			
			if (  is.null( x$usedObj[[uObj]])){
				x$usedObj[[uObj]] <- randomForest( x= t(as.matrix(fit_4_rf(x)$dat)), y=factor(x$samples[, group]),ntree=2000 )
			}
			gc()
			x$usedObj[[uObj]]
		}
)

