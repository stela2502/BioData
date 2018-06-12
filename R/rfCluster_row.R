#' @name rfCluster_row
#' @aliases 'rfCluster_row,BioData-method
#' @title rfCluster_row
#' @name rfCluster_row-methods
#' @docType methods
#' @description This fucntion uses the RFclust.SGE to create fandomForest based unsupervised clusters on a subset of the data.
#' @description Default is on 200 cells using all (provided) genes with 500 forests and 500 trees per forest for 5 repetitions.
#' @description You are asked to give a k numer of expected clusters (better too many than too little), classifies the total 
#' @description data using the 5 different unsupervised runs and all cluster ids from these runs are merged into the final cluster id.
#' @description This <summaryCol> will be part of the return objects samples table, together with a <usefulCol> where
#' @description all clusters with less than 10 cells have been merged into the 'gr. 0'.
#' @description The final results will be reported as new columns in the samples table containing the 'name'
#' @param x the single cells ngs object
#' @param email your email to use together with the SGE option
#' @param SGE whether to use the sun grid engine to calculate the rf grouping
#' @param rep how many repetitions for the random forest grouping should be run (default = 5)
#' @param slice how many processes should be started for each random forest clustering (default = 30)
#' @param bestColname the column name to store the results in
#' @param k the numer of expected clusters (metter more than to view)
#' @param subset how many cells should be randomly selected for the unsupervised clustering (default = 200)
#' @param name if you want to run multiple RFclusterings on e.g. using different input genes you need to specify a name (default ='RFclust')
#' @param nforest the numer of forests to grow for each rep (defualt = 500)
#' @param ntree the numer of trees per forest (default = 500)
#' @param settings slurm settings list(A, t and p) which allow to run the rf clustering on a slurm backend
#' @return a SingleCellsNGS object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
if ( ! isGeneric('rfCluster_row') ){ setGeneric('rfCluster_row',
		function ( x, rep=1, SGE=F, email="none", k=16, slice=4, subset=200,nforest=500, ntree=500, name='RFclust', settings=list()){
			standardGeneric('rfCluster_row')
		}
)
}else {
	print ("Onload warn generic function 'rfCluster_row' already defined - no overloading here!")
}
setMethod('rfCluster_row', signature = c ('BioData'),
		definition = function ( x, rep=1, SGE=F, email="none", k=16, slice=4, subset=200 ,nforest=500, ntree=1000, name='RFclust_row', settings=list()) {
			x$name <- str_replace_all( x$name, '\\s+', '_')
			summaryCol=paste( 'All_groups', name,sep='_')
			usefulCol=paste ('Usefull_groups',name, sep='_')
			n= paste(x$name, name,sep='_')
			m <- max(k)
			OPATH <- file.path( x$outpath,str_replace( x$name, '\\s', '_'))
			opath = file.path( OPATH,name,"RFclust.mp" )
			
			if ( ! dir.exists(OPATH)){
				dir.create( OPATH )
			}
			if ( ! dir.exists(file.path(OPATH, name )) ){
				dir.create(file.path(OPATH, name ) )
			}
			if ( ! dir.exists(file.path(OPATH, name, "RFclust.mp")) ){
				dir.create(file.path(OPATH, name,"RFclust.mp" ) )
			}
			processed = FALSE
			single_res_row <- paste('RFgrouping',name)
			for ( i in 1:rep) {
				tname = paste(n,i,sep='_')
				if ( is.null(x$usedObj[['rfExpressionSets_row']][[tname]]) ){
					i <- length(x$usedObj$rfObj_row)+i
					## start the calculations!
					if ( dir.exists(opath)){
						if ( opath == '' ) {
							stop( "Are you mad? Not giving me an tmp path to delete?")
						}
						system( paste('rm -f ',opath,"/*",tname,'*', sep='') )
					}else {
						dir.create( opath )
					}
					total <- nrow(x$dat)
					if ( total-subset <= 20  && rep > 1) {
						stop( paste( 'You have only', total, 'samples in this dataset and request to draw random',subset, "samples, which leaves less than 20 cells to draw on random!") )
					}
					else if ( total < subset ){
						stop ( paste("You can not ask for more than the max of",total, "samples in the test dataset!") )			
					}
					if ( is.null(x$usedObj[['rfExpressionSets_row']])){
						x$usedObj[['rfExpressionSets_row']] <- list()
						x$usedObj[['rfObj_row']] <- list()
					}
					if ( length( x$usedObj[['rfExpressionSets_row']] ) < i  ) {
						x$usedObj[['rfExpressionSets_row']][[ i ]] <- transpose(reduceTo( x,'row',to= rownames(x$dat)[sample(c(1:total),subset)], name=tname, copy=TRUE ))
						fit_4_rf(x$usedObj[['rfExpressionSets_row']][[ i ]], copy=F)
						if ( length(settings) > 0 ) {
							x$usedObj[['rfObj_row']][[ i ]] <- RFclust.SGE::RFclust.SGE ( 
									dat=as.data.frame(x$usedObj[['rfExpressionSets_row']][[ i ]]$data()), 
									SGE=F, slices=slice, email=email, tmp.path=opath, 
									name= tname, slurm=T,settings=settings 
							)
						}else {
							x$usedObj[['rfObj_row']][[ i ]] <- RFclust.SGE::RFclust.SGE ( 
									dat=as.data.frame(x$usedObj[['rfExpressionSets_row']][[ i ]]$data()), 
									SGE=SGE, slices=slice, email=email, tmp.path=opath, name= tname 
							)
						}
					}
					names(x$usedObj[['rfExpressionSets_row']])[i] <- tname
					names(x$usedObj[['rfObj_row']])[i] <- tname
					x$usedObj[['rfObj_row']][[ i ]] <- RFclust.SGE::runRFclust ( 
							x$usedObj[['rfObj_row']][[ i ]] , 
							nforest=nforest,
							ntree=ntree,
							name=tname 
					)
					if ( SGE){
						print ( "You should wait some time now to let the calculation finish! check: system('qstat -f') -> re-run the function")
					}
					else {
						print ( "You should wait some time now to let the calculation finish! -> re-run the function")
						print ( "check: system( 'ps -Af | grep Rcmd | grep -v grep')")
					}
				}
				else {
					i <- match( tname,names(x$usedObj$rfObj_row) )
					## read in the results
					try ( x$usedObj[['rfObj_row']][[ i ]] <- RFclust.SGE::runRFclust ( x$usedObj[['rfObj_row']][[ i]] , nforest=nforest, ntree=ntree, name=tname ) )
					if ( is.null(x$usedObj[["rfExpressionSets_row"]][[i]]$usedObj[['transposed']])){
						transpose( x$usedObj[["rfExpressionSets_row"]][[i]] ) ## changes the R6 object
					}
					if ( ! is.null(x$usedObj[['rfObj_row']][[ i ]]@RFfiles[[tname]]) ){
						stop( "please re-run this function later - the clustring process has not finished!")
					}
					for ( a in k ){
						# remember the x$usedObj[["rfExpressionSets_row"]][[i]] object is transposed!
						createRFgrouping_row( x, RFname=tname,  k=a, single_res_row = paste( single_res_row, a) )
						print ( paste("Done with cluster n=",a))
					}
										
					print ( paste("Done with RF clustering"))
					processed = TRUE
				}
			}
			invisible(x)		
		}
)




#' @name createRFgrouping_row
#' @aliases createRFgrouping_row,BioData-method
#' @rdname createRFgrouping_row-methods
#' @docType methods
#' @description Create a sample grouping data from one RFclust.SGE object
#' @param x the BioData object
#' @param RFname the name of the RFclust.SGE object in the BioData object. This object has to be populized with data!
#' @param k the number of wanted groups ( default = 10)
#' @param single_res_row the new column in the samples table default= paste('RFgrouping', RFname)
#' @param colFunc a function giving the colours back for the grouping (gets the amount of groups) default = function(x){rainbow(x)}
#' @title description of function createRFgrouping_row
#' @export 
if ( ! isGeneric('createRFgrouping_row') ){ setGeneric('createRFgrouping_row', ## Name
		function ( x, RFname='notExisting', k=10, single_res_row = paste('RFgrouping',RFname), colFunc=NULL) { ## Argumente der generischen Funktion
			standardGeneric('createRFgrouping_row') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
}else {
	print ("Onload warn generic function 'createRFgrouping_row' already defined - no overloading here!")
}

setMethod('createRFgrouping_row', signature = c ('BioData'),
		definition = function ( x, RFname='notExisting', k=10, single_res_row = paste('RFgrouping',RFname), colFunc=NULL) {
			if ( is.na( match( RFname, names(x$usedObj[['rfObj_row']])))){
				stop( paste("the RFname",RFname,"is not defined in this object; defined grouings are:",paste(names(x$usedObj[['rfObj_row']]), collapse=" ",sep=', ') ) )
			}
			groups <- createGroups( x$usedObj[['rfObj_row']][[RFname]], k=k, name=RFname )
			x$usedObj[['rfExpressionSets_row']][[RFname]]$samples <- 
					cbind ( x$usedObj[['rfExpressionSets_row']][[RFname]]$samples, groups[,3:(2+length(k))] )
			le <- ncol(x$usedObj[['rfExpressionSets_row']][[RFname]]$samples)
			colnames(x$usedObj[['rfExpressionSets_row']][[RFname]]$samples)[(le-length(k)+1):le] <- 
					paste('group n=',k)
			m <- max(k)
			## create the predictive random forest object
			if ( all.equal( sort(colnames(x$usedObj[['rfObj_row']][[RFname]]@dat)), sort(rownames(x$dat)) ) == TRUE ) {
				## use the column in grouping
				print ( "using the calcualted grouping")
				for ( id in 1:length(k) ){
					mat <- match(rownames(x$dat), colnames(x$usedObj[['rfObj_row']][[RFname]]@dat))
					x$annotation[, paste( single_res_row, ' n=', k[id], sep="") ] = factor(groups[mat,2+id], levels=c(1:k[id]))
					x <- colors_4( x, paste( single_res_row, ' n=', k[id], sep="")  )	
				}
			}
			else {
				#predict based on the RFdata
				print ( "predicting on the calculated grouping" )
				
				x$usedObj[['rfExpressionSets_row']][[RFname]] <- 
						bestGrouping( x$usedObj[['rfExpressionSets_row']][[RFname]], 
								group=paste('group n=', m), 
								bestColname = paste('OptimalGrouping',m ,RFname)
						)
				print( "rf predict")
				
				x$annotation[, paste( single_res_row) ] <-
						predict( 
								x$usedObj[['rfExpressionSets_row']][[RFname]]$usedObj[[paste( 'predictive RFobj group n=',m) ]], 
								as.matrix(fit_4_rf(x)$dat)
						)
				x$annotation[, paste( single_res_row) ] <- factor( x$annotation[, paste( single_res_row) ], levels= 1:m )
				x <- colors_4( x, single_res_row )
			}
			gc()
			invisible(x)
		} 
)


