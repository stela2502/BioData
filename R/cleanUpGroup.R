#' @name cleanUpGroup
#' @aliases cleanUpGroup,BioData-method
#' @rdname cleanUpGroup-methods
#' @docType methods
#' @description Clean up intermixed groups
#' @param x the BioData object
#' @param group the samples column to clean up
#' @param otherGroup the samples column to look for contaminations
#' @param max_cells how many cells to max use to create the predictive RF object (default = 10 [cells])
#' @param min_cells if two groups inside one cluster exceed min_cells one more cell group is created (defualt =1\% of all cells)
#' @param ... additional variables for the randomForest call
#' @title description of function cleanUpGroup
#' @export 
setGeneric('cleanUpGroup', ## Name
		function ( x, group, otherGroup, max_cells=10, min_cells=10, ...) { ## Argumente der generischen Funktion
			standardGeneric('cleanUpGroup') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('cleanUpGroup', signature = c ('BioData'),
		definition = function ( x, group, otherGroup, max_cells=10, min_cells=NULL, ... ) {
			if ( is.null(min_cells)) { min_cells=round( ncol(x$dat) / 100) }
			## get the main otherGroup for each group 
			m=matrix(0, nrow=2, ncol=2)
			rownames(m) = c("g2", "g1")
			colnames(m) = c("S1", "S2")
			h = as_BioData( m )
			r = lapply( levels(x$samples[, group]),
					function(n, h ) {
						t = table(x$samples[, otherGroup][which(x$samples[, group] == n )])
						Groups= NULL
						if ( length(which( t > min_cells)) > 0 ) {
							Groups = names(t)[which(t > min_cells)]
						}else {
							Groups = names(t[which( t == max(t))])[1]
						}
						
						for ( mainOther in Groups ){
							groupIDs <- which(x$samples[, group] == n )
							OK = groupIDs[ which(x$samples[ groupIDs ,otherGroup ] ==  mainOther )]
							if ( length(OK) >= max_cells) {
								OK = sample(OK, max_cells)
							}
							#OK = paste(collapse=";", OK )
							if ( is.null(h$usedObj$save) ){
								h$usedObj$save= as.matrix( data.frame( 
												rep(mainOther,length(OK)), 
												rep(paste( n, mainOther, sep="_"),
														length(OK)
											), OK) )
							}else {
								h$usedObj$save = rbind(h$usedObj$save, as.matrix( data.frame( rep(mainOther,length(OK)), 
														rep(paste( n, mainOther, sep="_"),length(OK)), OK) ) )
							}
						}
						NULL
					} , h
			)
		
			r= h$usedObj$save
			rownames(r) = NULL
			r = cbind(r, x$samples[as.numeric(r[,3]), x$sampleNamesCol])
			OK = reduceTo(x, what='col', to=colnames(x$dat)[as.numeric(r[,3])], name="OK" , copy=T)
			OK$samples$New_Grouping = as.factor(r[,2])
			print ("Create predictive RF object")
			#browser()
			RFobj = bestGrouping(OK, 'New_Grouping', ...)
			x$samples[,paste(group, 'original')] =  x$samples[,group]
			x$samples[, group] = factor( predict( RFobj , t(as.matrix(x$data())) ) )#, levels=levels(x$samples[, group]))
			group_intersect_order( x,group, paste(group, 'original') )
			colors_4( x, group, force=T)
			print( paste("Grouping", group,"cleaned by grouping", otherGroup ))
			invisible(x)
		}  )
