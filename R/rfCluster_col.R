#' @name rfCluster_col
#' @aliases 'rfCluster_col,BioData-method
#' @title rfCluster_col
#' @name rfCluster_col-methods
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
#' @param rep how many repetitions for the random forest grouping should be run (default = 1)
#' @param slice how many processes should be started for each random forest clustering (default = 4)
#' @param k the numer of expected clusters (better more than to view)
#' @param subset how many cells should be randomly selected for the unsupervised clustering (default = 200)
#' @param name if you want to run multiple RFclusterings on e.g. using different input genes you need to specify a name (default ='RFclust')
#' @param nforest the numer of forests to grow for each rep (defualt = 500)
#' @param ntree the numer of trees per forest (default = 500)
#' @param settings slurm settings list(A, t and p) which allow to run the rf clustering on a slurm backend
#' @param ids the ids for a subset of samples to be analyzed (default NULL)
#' @return a SingleCellsNGS object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
if ( ! isGeneric('rfCluster_col') ){ methods::setGeneric('rfCluster_col',
			function ( x, rep=1, SGE=F, email='none', k=16, slice=4, 
					subset=200,nforest=500, ntree=500, name='RFclust',
					settings=list(), ids=NULL){
				standardGeneric('rfCluster_col')
			}
	)
}else {
	print ("Onload warn generic function 'rfCluster_col' already defined - no overloading here!")
}
setMethod('rfCluster_col', signature = c ('BioData'),
		definition = function ( x, rep=1, SGE=F, email="none", k=16, slice=4, 
				subset=200 ,nforest=500, ntree=1000, name='RFclust',
				settings=list(), ids=NULL) {
			
			if ( rep > 1) {
				lapply(1:rep, function(i) { 
							rfCluster_col(
									x, rep=1, SGE=SGE, email=email, 
									k=k, slice=slice, subset=subset ,nforest=nforest, 
									ntree=ntree, name=paste(sep="_",name, i ), settings=settings, ids= ids ) 
						} )
			}
			x$name <- stringr::str_replace_all( x$name, '\\s+', '_')
			summaryCol=paste( 'All_groups', name,sep='_')
			usefulCol=paste ('Usefull_groups',name, sep='_')
			n= paste(x$name, name,sep='_')
			m <- max(k)
			OPATH <- file.path( x$outpath,stringr::str_replace( x$name, '\\s', '_'))
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
			single_res_col <- paste('RFgrouping',name)
			i = 1
			tname = paste(n,i,sep='_')
			if ( is.null(x$usedObj[['rfExpressionSets']][[tname]]) ){
				i <- length(x$usedObj$rfObj)+i
				## start the calculations!
				if ( dir.exists(opath)){
					if ( opath == '' ) {
						stop( "For security reasons you must provide an absolute tmp path - files in there will be deleted?")
					}
					system( paste('rm -f ',opath,"/*",tname,'*', sep='') )
				}else {
					dir.create( opath )
				}
				total <- ncol(x$dat)
				if ( total-subset <= 20  && rep > 1) {
					stop( paste( 'You have only', total, 'samples in this dataset and request to draw random',subset, "samples, which leaves less than 20 cells to draw on random!") )
				}
				else if ( total < subset ){
					stop ( paste("You can not ask for more than the max of",total, "samples in the test dataset!") )			
				}
				if ( is.null(x$usedObj[['rfExpressionSets']])){
					x$usedObj[['rfExpressionSets']] <- list()
					x$usedObj[['rfObj']][[ i ]] <- list()
				}
				
				if ( length( x$usedObj[['rfExpressionSets']] ) < i  ) {
					if ( ! is.null(ids) ) {
						x$usedObj[['rfExpressionSets']][[ i ]] <- 
						  reduceTo( x, what='col', to=colnames(x$dat)[ids], name=tname, copy=TRUE )
					}else {
						x$usedObj[['rfExpressionSets']][[ i ]] <- 
						  reduceTo( x, what='col', to=colnames(x$dat)[sample(c(1:total),subset)], name=tname, copy=TRUE )
					}
					## here I need to get rid of the -1 values!
					fit_4_rf(x$usedObj[['rfExpressionSets']][[ i ]], copy=F)
					if ( length(settings) > 0 ) {
						#browser()
						x$usedObj[['rfObj']][[ i ]] <- RFclust.SGE::RFclust.SGE ( 
								dat=x$usedObj[['rfExpressionSets']][[ i ]]$data(), 
								SGE=F, slices=slice, email=email, tmp.path=opath, 
								name= tname, settings=settings, slurm=T 
						)
					}else {
						x$usedObj[['rfObj']][[ i ]] <- RFclust.SGE::RFclust.SGE ( 
								dat=x$usedObj[['rfExpressionSets']][[ i ]]$data(), 
								SGE=SGE, slices=slice, email=email, tmp.path=opath, name= tname 
						)
					}
				}
				names(x$usedObj[['rfExpressionSets']]) [i] <- tname
				names(x$usedObj[['rfObj']]) [i] <- tname
				x$usedObj[['rfObj']][[ i ]] <- RFclust.SGE::runRFclust ( x$usedObj[['rfObj']][[ i ]] , nforest=nforest, ntree=ntree, name=tname )
				if ( SGE){
					print ( "You should wait some time now to let the calculation finish! check: system('qstat -f') -> re-run the function")
				}
				else {
					print ( "You should wait some time now to let the calculation finish! -> re-run the function")
					print ( "check: system( 'ps -Af | grep \"R.*BATCH\" | grep -v grep')")
				}
			}
			else {
				i <- match( tname,names(x$usedObj$rfObj) )
				## read in the results
				try ( x$usedObj[['rfObj']][[ i ]] <- runRFclust ( x$usedObj[['rfObj']][[ i]] , nforest=nforest, ntree=ntree, name=tname ) )
				if ( ! is.null(x$usedObj[['rfObj']][[ i ]]@RFfiles[[tname]]) ){
					stop( "please re-run this function later - the clustring process has not finished!")
				}
				for ( a in k ){
					x$usedObj[["rfExpressionSets"]][[i]]$samples <- 
							x$usedObj[["rfExpressionSets"]][[i]]$samples[ ,
									is.na(match ( colnames(x$usedObj[["rfExpressionSets"]][[i]]$samples), paste('group n=',a) ))==T 
							]
				}
				x <- createRFgrouping_col( x, RFname=tname,  k=k, single_res_col = paste( single_res_col, i) )
				
				print ( paste("Done with cluster",i))
				processed = TRUE
			}
			
			gc()
			
			invisible(x)		
		}
)




#' @name createRFgrouping_col
#' @aliases createRFgrouping_col,BioData-method
#' @rdname createRFgrouping_col-methods
#' @docType methods
#' @description Create a sample grouping data from one RFclust.SGE object
#' @param x the BioData object
#' @param RFname the name of the RFclust.SGE object in the BioData object. This object has to be populized with data!
#' @param k the number of wanted groups ( default = 10)
#' @param single_res_col the new column in the samples table default= paste('RFgrouping', RFname)
#' @param colFunc a function giving the colours back for the grouping (gets the amount of groups) default = function(x){rainbow(x)}
#' @title description of function createRFgrouping_col
#' @export 
if ( ! isGeneric('createRFgrouping_col') ){ methods::setGeneric('createRFgrouping_col', ## Name
			function ( x, RFname, k=10, single_res_col = paste('BioData',RFname), colFunc=NULL) { 
				standardGeneric('createRFgrouping_col')
			}
	)
}else {
	print ("Onload warn generic function 'createRFgrouping_col' already defined - no overloading here!")
}


setMethod('createRFgrouping_col', signature = c ('BioData'),
		definition = function ( x, RFname, k=10, single_res_col = paste('RFgrouping',RFname), colFunc=NULL) {
			if ( is.na( match( RFname, names(x$usedObj[['rfObj']])))){
				stop( paste("the RFname",RFname,"is not defined in this object; defined grouings are:",paste(names(x$usedObj[['rfObj']]), collapse=" ",sep=', ') ) )
			}
			groups <- RFclust.SGE::createGroups( x$usedObj[['rfObj']][[RFname]], k=k, name=RFname )
			## store the MDS representation of the rfClust dissimilarity object
			## in case all samples have been used to create the object.
			if ( is.null (x$usedObj$MDS[[RFname]]) & nrow(x$usedObj$rfObj[[RFname]]@distRF[[RFname]]) == ncol(x$dat) ) {
				a <- (x$usedObj$rfObj[[RFname]]@distRF[[RFname]])
				d <- stats::cmdscale(a,3)
				m <- match( colnames(x$dat), rownames(d) )
				x$usedObj$MDS[[single_res_col]] <- d[m,]
			}
			x$usedObj[['rfExpressionSets']][[RFname]]$samples <- 
					cbind ( x$usedObj[['rfExpressionSets']][[RFname]]$samples, groups[,3:(2+length(k))] )
			
			le <- ncol(x$usedObj[['rfExpressionSets']][[RFname]]$samples)
			colnames(x$usedObj[['rfExpressionSets']][[RFname]]$samples)[(le-length(k)+1):le] <- 
					paste('group n=',k)
			m <- max(k)
			## create the predictive random forest object
			browser()
			if ( all.equal(sort( colnames(x$usedObj[['rfObj']][[RFname]]@dat) ), sort( colnames(x$dat) ) ) == TRUE ) {
				## use the column in grouping
				print ( "using the calcualted grouping")
				mat <- match(colnames(x$dat), colnames(x$usedObj[['rfObj']][[RFname]]@dat))
				for ( id in 1:length(k) ){
					#browser()
					x$samples[, paste( single_res_col, ' n=', k[id], sep="") ] = factor(groups[mat,2+id], levels=c(1:k[id]))
					lapply( paste( single_res_col, ' n=', k[id], sep=""), function( n ) { colors_4(x,n) }) 
				}
			}else {
				#predict based on the RFdata
				print ( "predicting on the calculated grouping" )
				RFobj <- bestGrouping( x$usedObj[['rfExpressionSets']][[RFname]], group=paste('group n=', m), bestColname = paste('OptimalGrouping',m ,RFname))
				print( "rf predict")
				x$samples[, paste( single_res_col) ] <-	stats::predict( RFobj , t(as.matrix(x$data())) )
				x$samples[, paste( single_res_col) ] <- factor( x$samples[, paste( single_res_col) ], levels= 1:m )
				x <- colors_4( x, single_res_col )
			}
			invisible(x)
		} 
)


