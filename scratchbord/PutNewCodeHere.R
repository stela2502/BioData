combineGroupings <- function( x, colGroups, new_name="Merged Group 1", minCellOverlap=5 , minGroupOverlap=length(colGroups)-2, minCellsInReturnGroup=NULL ) {
	
	info = as.matrix( x$samples[,colGroups] )
	
	
	
	if ( is.null(minCellsInReturnGroup)){
		minCellsInReturnGroup =
			max(
				unlist(lapply( levels(merged$samples$sname), 
					function(n) { 
						length( which( merged$samples$sname == n)) 
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
	
	print ( paste("remapping", length(IDS),"small group cells with lowered group overlap to ", minGroupOverlap) )
	
	steps = ceiling(length(IDS)/100)
	pb <- progress_estimated(100)
	a= 0;
	for ( i in IDS ) {
		closestCells(i, orNULL=TRUE)
		a = a+1
		if ( a %% steps == 0 ) {
			pb$tick()$print()
			#print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
		}
	}
	pb$stop()
	#pb$tick()$print()
	
	x$usedObj$tmpres[which(is.na(x$usedObj$tmpres))] = 0
	
	x$usedObj$tmpgroupID = max(as.numeric(x$usedObj$tmpres))
	rname = which( as.numeric(x$usedObj$tmpres) == 0 )
	print ( paste( length(ungroup), "tiny groups with less than" , 
					minCellsInReturnGroup,"cells per group have been combined into the 'ungrouped' group (n=",
					length(rname),")"
			))
	if ( length(rname) > 0 ){
		x$usedObj$tmpres[rname] = "ungrouped"
		rname = c(1:x$usedObj$tmpgroupID, "ungrouped")
	}else {
		rname = c(1:x$usedObj$tmpgroupID)
	}
	
	x$samples[,new_name] = factor(x$usedObj$tmpres, levels=rname)
	x$usedObj$tmpres = NULL
	x$usedObj$tmpgroupID = NULL
	res
}
