#' @name combineGroupings
#' @aliases combineGroupings,BioData-method
#' @rdname combineGroupings-methods
#' @docType methods
#' @description combine different grouping informations to create overlaps. 
#' Implemented to combine different subset RFclusterings into one including all features.
#' @param x The BioData object
#' @param colGroups which groupings to combine (samples columns)
#' @param new_name the name of the new grouping (default="Merged Group 1" will replace existing information)
#' @param minCellOverlap while calculating the grouing overlap - how many cells need to overlap with each other default=5 
#' @param minGroupOverlap while caulculating the overlap - how many groups should overlap in order to considder a cell overlapping default=length(colGroups)-2
#' @param minCellsInReturnGroup merge all groups with less than this number of cells into a ungrouped group default=1\% of the bigges sname subgroup
#' @title description of function combineGroupings
#' @details This function will calculate the overlap in group names for each cell to all other cells,
#' @details select the cells that match the requirement ( minCellOverlap and minGroupOverlap) and asigne the group id of the biggest group
#' @details to this cell.
#' @details The same algorythm will be applied for all cells that end up in small groups before the final ungrouped cells are assigned.
#' @examples combineGroupings(BioDataObj, colnames(BioDataObject$samples)[1:10], minCellOverlap = 5, minGroupOverlap = 8, minCellsInReturnGroup = 20)
#' @export 
setGeneric('combineGroupings', ## Name
	function ( x, colGroups, new_name="Merged Group 1", minCellOverlap=5 , minGroupOverlap=length(colGroups)-2, minCellsInReturnGroup=NULL ) { 
		standardGeneric('combineGroupings')
	}
)

setMethod('combineGroupings', signature = c ('BioData'),
	definition = function ( x, colGroups, new_name="Merged Group 1", minCellOverlap=5 , minGroupOverlap=length(colGroups)-2, minCellsInReturnGroup=NULL ) {
	
	info = as.matrix( x$samples[,colGroups] )
	
	
	
	if ( is.null(minCellsInReturnGroup)){
		minCellsInReturnGroup =
			max(
				unlist(lapply( levels(x$samples$sname), 
					function(n) { 
						length( which( x$samples$sname == n)) 
					} )) / 100
		)
		print (paste("minCellsInReturnGroup set to",minCellsInReturnGroup) )
	}
	
	closestCells = function(colId, orNULL=FALSE) {
		my = info[colId,]
		d = apply (info[-colId,] ,1, function (v) {
					sum ( my == v )
				} )
		names(d) = colnames(x$dat)[-colId]

		t <- table(d) ## number of cells with names(t) overlap
		## I want to get the number of overlapping groups per cell with the most cells from t
		m = NULL
		if ( length( which(as.numeric(names(t)) >= minGroupOverlap & t > minCellOverlap ) ) > 0 ){
			## take the most close cells only
			useful = t[which(as.numeric(names(t)) >= minGroupOverlap & t > minCellOverlap )]
			m = as.numeric(names(useful)[length(useful)])
		}
		else {
			m = as.numeric(names(t)[length(t)]) ## get the highest numer of overlaps
			if (m < minGroupOverlap) {
				m = minGroupOverlap
			}
		}
		x$usedObj$tmpM[i] = m
		r = names(which( d>= m ))
		ok = sort(match( sort(r), colnames(x$dat)))
		
		if (any(is.na(x$usedObj$tmpres[ok]) == FALSE) ){
			t = table( x$usedObj$tmpres[ok] )
			## set the Id to the one with the most entries!
			
			highest_id = order(t,decreasing=T)[1]
			x$usedObj$tmpN[i] = t[highest_id] ## store the number of cells supporting this group
			x$usedObj$tmpres[i] = as.numeric(names(t)[highest_id])
			#	print ( paste( i, "set to", x$usedObj$tmpres[i]) )
		}else {
			if ( orNULL ) {
				x$usedObj$tmpres[i] = 0
			}else {
				x$usedObj$tmpgroupID=x$usedObj$tmpgroupID+1
				x$usedObj$tmpres[i] = x$usedObj$tmpgroupID
			}
		}
		NULL
		#ok ## return the col ids of the best match
	}
	
	if ( is.null(x$usedObj$tmpres)){
		
		x$usedObj$tmpres =  rep( NA , ncol(x$dat) )
		x$usedObj$tmpM =  rep( NA , ncol(x$dat) )
		x$usedObj$tmpN =  rep( NA , ncol(x$dat) )
		
		x$usedObj$tmpgroupID=0
		steps = ceiling(ncol(x$dat)/100)
		pb <- progress_estimated(100)
		print ( paste( "Calculating group overlap for", ncol(x$dat),"cells" ) )
		for ( i in 1:ncol(x$dat) ) {
			closestCells(i)
			if ( i %% steps == 0 ) {
				pb$tick()$print()
				#print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
			}
		}
		#pb$tick()$print()
		pb$stop()
	}
	res = x$usedObj$tmpres
	
	print ( paste( "created",x$usedObj$tmpgroupID, "temporary groups" ))
	
	#clear out tiny groups
	t <- table(as.numeric(x$usedObj$tmpres))
	ungroup = names(which ( t < minCellsInReturnGroup ))
	#x$usedObj$tmpres = as.numeric(x$usedObj$tmpres)
	for ( i in rev(as.numeric(ungroup)) ) {
		reduce = which(as.numeric(x$usedObj$tmpres) > i )
		## has it been a position in the order that keeps the cell from finding a good group?
		IDS = which( as.numeric(x$usedObj$tmpres)==i )
		if ( length(IDS) > 0 ){
			#browser()
			x$usedObj$tmpres[IDS] = 0 ## this group has to be erased
		}
		x$usedObj$tmpres[reduce] = x$usedObj$tmpres[reduce] -1
	}
	
	## So now we have a lot of 0 groups and I would like to re-map them
	x$usedObj$tmpgroupID=max(as.numeric(x$usedObj$tmpres))
	IDS = which( x$usedObj$tmpres == 0 )
	minGroupOverlap = minGroupOverlap -2
	
	## useless - use the others to build the module to include them all.
	OK = reduceTo( x, copy=T, to=colnames(x$dat)[which( x$usedObj$tmpres != 0)], what='col' )
	OK$samples$tmp_Group = x$usedObj$tmpres[which( x$usedObj$tmpres != 0)]
	
	print ( paste("remapping", length(IDS),"small group cells with lowered group overlap to ", minGroupOverlap) )
	print("Creating predictive RFobject")
	RFobj = bestGrouping( OK, 'tmp_Group')
	print ("apply grouping on the whole dataset")
	x$samples[, new_name] <- factor(predict( RFobj , t(as.matrix(x$data())) ) )
	colors_4(x, new_name)	
	print ( paste("remapping", length(IDS),"small group cells with lowered group overlap to ", minGroupOverlap) )

	res
} )
