#' @name mds
#' @aliases mds,BioData-method
#' @rdname mds-methods
#' @docType methods
#' @description Calculates the MDS for a given MDS type and stores the 3 dimensions 
#' in the object for later use. All projections use a n=100 PCA projection as starting material instead of the raw data.
#' @param dataObj the BioData object
#' @param mds.type Which MDS function should be called default="PCA"
#' @param onwhat condense which dataset at the moment only Expression is supported default='expression'
#' @param genes do it on genes not on samples (default = F)
#' @param LLEK the neighbours in the LLE algorithm (default=2)
#' @param useRaw base the projection on the raw data and not the n=100 PCA data (default FALSE)
#' @param pythonEnv if needed the virtual environment for the MulticoreTSNE call set to something like 'source <myPrivateEnvPath>/bin/activate'  (default = NULL)
#' @param dim the number of dimensions to return ( default 3)
#' @title Calculate MDS projections for the 3D Make3Dobj function
#' @export 
if ( ! isGeneric('mds') ){ setGeneric('mds', ## Name
			function ( dataObj, ..., mds.type="PCA" , onwhat ='Expression', genes=F,  LLEK=2, useRaw=F, pythonEnv=NULL, dim=3) { ## Argumente der generischen Funktion
				standardGeneric('mds') ## der Aufruf von standardGeneric sorgt für das Dispatching
			}
	)
}else {
	print ("Onload warn generic function 'mds' already defined - no overloading here!")
}

setMethod('mds', signature = c ('BioData'),
		definition = function ( dataObj, ..., mds.type="PCA", onwhat ='Expression', genes=F, LLEK=2, useRaw=F, pythonEnv=NULL, dim=3 ) {
			## the code is crap re-code!!
			mds_store = NULL
			if ( mds.type=="PCA" ) {
				seRaw=F
			}
			if ( onwhat == 'Expression'){
				n = ncol(dataObj$data())
				if ( n > 100) 
					n = 100
				rerun = 0
				if ( is.null(dataObj$usedObj$pr) ){
					rerun = 1	
				}else if ( all.equal( rownames(dataObj$usedObj$pr@scores), colnames(dataObj$dat) ) == F ) {
					rerun = 1
				}
				if ( rerun == 1) {
					tmp = dataObj$data()
					bad = which(tmp@x == -1)
					if ( length(bad) > 0 ) {
						tmp[bad] = 0
					}
					dataObj$usedObj$pr <- bpca( t(as.matrix(tmp)), nPcs=100 )
					
					#dataObj$usedObj$pr <- irlba::prcomp_irlba ( t(tmp), center=T, n=n )
					
					rm(tmp)
					#rownames(dataObj$usedObj$pr$x) = colnames(dataObj$dat)
					rownames(dataObj$usedObj$pr@scores) = colnames(dataObj$dat)
				}
			}else {
				stop( paste("Sorry, the option onwhat",onwhat,"is not supported") )
			}
			
			if ( useRaw ) {
				mds_store <- 'MDS'
				tab=t(as.matrix(dataObj$data()))
			}else {
				mds_store <- 'MDS_PCA100'
				#tab <- dataObj$usedObj$pr$x
				tab <- dataObj$usedObj$pr@scores
			}
			
			if ( genes ) {
				n = ncol(dataObj$data())
				if ( n > 100) 
					n = 100
				if ( is.null(dataObj$usedObj$prGenes) |  all.equal( rownames(dataObj$usedObj$prGenes@scores), rownames(dataObj$dat) ) == F ) {
					tmp = dataObj$data()
					bad = which(tmp@x  == -1)
					if ( length(bad) > 0 ) {
						tmp[bad] = 0
					}
					dataObj$usedObj$prGenes <- bpca( as.matrix(tmp), nPcs=100 )
					#dataObj$usedObj$prGenes <- irlba::prcomp_irlba ( tmp, center=T, n=n )
					rm(tmp)
					#rownames(dataObj$usedObj$prGenes$x) = rownames(dataObj$dat)
					rownames(dataObj$usedObj$prGenes@scores) = rownames(dataObj$dat)
				}
				if ( useRaw ) {
					mds_store <- 'MDSgene'
					tab <- t(tab)
				}else {
					mds_store <- 'MDSgene_PCA100'
					tab <- dataObj$usedObj$prGenes@scores	
				}	
			}
			if ( dim != 3){
				mds_store = paste( mds_store, 'dim',dim, sep="_")
			}
			
			## define mds storage position
			if ( is.null( dataObj$usedObj[[mds_store]]) ) {
				dataObj$usedObj[[mds_store]] = list()
			}
			this.k <- paste(onwhat,mds.type)
			
			## MDS code
			if ( (is.null(dataObj$usedObj[[mds_store]][[this.k]])) ||  all.equal( rownames(dataObj$usedObj[[mds_store]][[this.k]]), colnames(dataObj$dat) )==F ) {
				mds.proj <- NULL
				pr <- NULL
				#system ( 'rm loadings.png' )
				
				if(mds.type == "PCA"){
					mds.proj <- tab[,1:dim]
					try( {
								png ( file=file.path( dataObj$outpath,'loadings.png'), width=1000, height=1000 )
								plot (  pr$rotation[,1:2] , col='white' );
								text( pr$rotation[,1:2], labels= rownames(pr$rotation), cex=1.5 )
								dev.off()
								write.table( cbind( Genes = rownames(pr$rotation), pr$rotation[,1:2] ), 
										file=file.path( dataObj$outpath,'gene_loadings.xls') , row.names=F, sep='\t',quote=F )
							})
				}
				#	mds.trans <- prcomp(t(tab))$x[,1:3]
				
			} else if ( mds.type=='DM') {
				if (!library("destiny", quietly = TRUE,logical.return=TRUE )) {
					stop("package 'destiny' needed for this function to work. Please install it.",
							call. = FALSE)
				}
				if ( ! exists('sigma', mode='numeric') ){
					if ( exists('find_sigmas', where='package:destiny', mode='function') ){
						sigmas <- destiny::find_sigmas(tab, verbose=F)
						sigma <- destiny::optimal_sigma(sigmas)
					}else {
						sigmas <- destiny::find.sigmas(tab, verbose=F)
						sigma <- destiny::optimal.sigma(sigmas)
					}
					
				}
				if ( !exists('distance', mode='character')){
					distance = 'cosine'
				}
				dm <- destiny::DiffusionMap(tab, distance = distance, sigma = sigma, n_eigs = dim)
				mds.proj <- destiny::as.data.frame(dm)[,1:dim]
				
			}else if ( mds.type == "TSNE"){
				
				## there is an extremely efficient python implementation of that available.
				## lets try to use that:
				## first export the '\t' separated file
				if ( file.exists( 'runTSNE.py' ) & ! file.exists( 'tSNE_dim3_coods.csv') ){
					print ( "The pyton process should be started and processing the data - please check manually")
					return (invisible(dataObj))
					
				}else if ( ! file.exists( 'tSNE_dim3_coods.csv' ) )  {
					
					write.table( t(tab), sep=",", file="TSNE_data.csv", quote=F )
					fileConn<-file( 'runTSNE.py' )
					writeLines(c(
									"from MulticoreTSNE import MulticoreTSNE as TSNE",
									"import pandas as pd",
									"import numpy as np",
									"",
									paste( sep="","tsne = TSNE(n_jobs=4,n_components=",dim,")"),
									"X = pd.io.parsers.read_csv('TSNE_data.csv',sep=',',index_col=0)",
									"",
									"Y = tsne.fit_transform(np.transpose(X))",
									"",
									"df=pd.DataFrame(data=Y,index=None)",
									"df.to_csv('tSNE_dim3_coods.csv')"
							), fileConn )
					close(fileConn)
					if (! is.null(pythonEnv)) {
						system( paste( pythonEnv ,"&&", 'python' , 'runTSNE.py' ) )					
					}else {
						system( paste(Sys.which('python'), 'runTSNE.py' ) )
					}
					print ("the external python script has been run - rerun this function to check if it is finished.")
					print ( "in case the python script does not produce output (1) try to install MulticoreTSNE grom its git resource or (2) use mds.type='TSNE_R'")
					return (invisible(dataObj))
				} else {
					## OK the output file has been produced
					mds.proj <- read.delim( 'tSNE_dim3_coods.csv', sep="," );
					rownames(mds.proj) <- rownames(tab)
				}
			}else if ( mds.type == "TSNE_R"){
				if (!library("Rtsne", quietly = TRUE,logical.return=TRUE )) {
					stop("package 'Rtsne' needed for this function to work. Please install it.",
							call. = FALSE)
				}
				if ( useRaw ){
					mds.proj <- Rtsne( tab, dims=dim , check_duplicates =F, pca_center=F, verbose=T, pca=T )$Y
				}else {
					## The data is already PCAed
					mds.proj <- Rtsne( tab, dims=dim , check_duplicates =F,  verbose=T, pca=F )$Y
				}
				
				rownames(mds.proj) <- rownames(tab)
				
			}
			
			else if ( mds.type == "LLE"){
				mds.proj <- LLE( tab, dim = dim, k = as.numeric(LLEK) )
				#	mds.trans <- LLE( t(tab), dim = 3, k = as.numeric(LLEK) )
				
			}else if ( mds.type == "ISOMAP"){
				mds.proj <- Isomap( tab, dim = dim, k = as.numeric(LLEK) )$dim3
				#	mds.trans <- Isomap( t(tab), dim = 3, k = as.numeric(LLEK) )$dim3
				
			}else if ( mds.type == "ZIFA" ) {
				#stop( "Sorry ZIFA has to be double checked - are you working on normalized data - than ZIFA can not be applied!")
				print ( "Running external python script to apply ZIFA dimensional reduction (Expression data only)" )
				if ( genes ) {
					write.table( dataObj$dat, file="ZIFA_input.dat", sep=" ", col.names=F, row.names=F , quote=F)
				}else {
					write.table( t(dataObj$dat), file="ZIFA_input.dat", sep=" ", col.names=F, row.names=F , quote=F)
				}
				write( c("from ZIFA import ZIFA","from ZIFA import block_ZIFA", "import numpy as np",
								"Y = np.loadtxt('ZIFA_input.dat')", "Z, model_params = ZIFA.fitModel( Y, 3 )", 
								"np.savetxt('TheMDS_ZIFA.xls', Z )" ), 
						file= 'ZIFA_calc.py' )
				system( "python ZIFA_calc.py" )
				Sys.sleep(5)
				mds.proj <- read.delim( "TheMDS_ZIFA.xls", sep=' ', header=F)
				if ( genes ) {
					rownames(mds.proj) <- rownames(dataObj$dat)
				}else {
					rownames(mds.proj) <- colnames(dataObj$dat)
				}
				colnames(mds.proj) <- c( 'x','y','z')
				
			} else if ( mds.type == "DDRTree" ) {
				DDRTree_res <- DDRTree( t(tab), dimensions=3)
				mds.proj <- t(DDRTree_res$Z)
				rownames(mds.proj) <- rownames(tab)
				dataObj$usedObj$DRRTree.genes <- DDRTree_res
				
			}
			else {
				print( paste("Sory I can not work on the mds.type option",mds.type) )
			}
			if ( genes ) {
				rownames(mds.proj) <- rownames(dataObj$dat)
			}else{
				rownames(mds.proj) <- make.names(colnames(dataObj$dat))
			}
			
			dataObj$usedObj[[mds_store]][[this.k]]<- mds.proj
			
			
			invisible(dataObj)
		} 
		)
