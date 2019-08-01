#' @name clusters
#' @aliases clusters,BioData-method
#' @rdname clusters-methods
#' @docType methods
#' @description This function clusters the data based on hclust, mclust or the kmeans algorithm.
#' This function can be used to cluster any data slot the BioData class offers (see combination of onwhat and clusterby).
#' @param dataObj the BioData object
#' @param clusterby is used to specify the MDS type if onwhat is set to MDS (PCA, TSNE_R ...) or the raw, zscored or raw slot for expression data (default 'raw')
#' @param useGrouping do nothing and simply use this grouping default=NULL
#' @param groups.n how many groups should be detected default= 3
#' @param ctype cluster type - either 'hierarchical clust', 'mclust' or 'kmeans' default = 'hierarchical clust'
#' @param onwhat This option selects the source data ('Expression', or 'MDS')
#' @param cmethod the method to used with the hclust clustering (default = 'ward.D2')
#' @param name the name for the new grouping (default = 'auto_clusters.1:n')
#' @title create a grouping based on either the raw data or a MDS projection
#' @examples
#' clusters( TestData , clusterby = "TSNE_R", groups.n = 3, ctype = "kmeans", onwhat= 'MDS', name = "kmeansTSNE_R clusters" )
#' @export

# if ( ! isGeneric('clusters') ){ 
setGeneric('clusters', ## Name
	function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3, 
			ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D2', name=NULL ) {
		standardGeneric('clusters')
	}
)

#}else {
#	print ("Onload warn generic function 'clusters' already defined - no overloading here!")
#}

setMethod('clusters', signature = c ('BioData'),
	definition = function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3,
				ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D2', name=NULL ) {
	
			clusters <- NULL
			hc <- NULL
			MDS_TYPE <- NULL
			if(onwhat=="Expression" & clusterby=="raw"){
				tab <- as.matrix(dataObj$data())
				m <- min(tab)
				if ( m == -1 | m == -21 ){ ## get rid of them
					bad = which( tab == m)
					if ( length(bad) > 0 ){
						tab[bad] = m +1
					}
				}
			}else if ( onwhat=="MDS" ) {
				check_and_replace <- function( name, list) {
					if ( length(grep(name, names(list) )) ==1 ){
						name = names(list)[grep(name, names(list) )]
						print ( paste( "name changed to",name))
					}
					name
				}
				MDS_TYPE= 'MDS'
				test = check_and_replace(clusterby, dataObj$usedObj[[MDS_TYPE]] )
				if ( is.na( match( test, names(dataObj$usedObj[[MDS_TYPE]] )))){
					MDS_TYPE= 'MDS_PCA100'
					test = check_and_replace(clusterby, dataObj$usedObj[[MDS_TYPE]] )
				}
				cn <- names(dataObj$usedObj[[MDS_TYPE]])
				if ( length(grep(test, cn  )) == 0) {
					dataObj <- mds( dataObj, onwhat="Expression", mds.type=clusterby)
				}else if (length(grep(test, cn  )) > 1 ) {
					stop( paste("Sorry, but I have more than one MDS type like = '",clusterby,"' :",paste(sep=",",cn[grep(clusterby, cn  )] ), sep="") )
				}
				clusterby = check_and_replace(clusterby, dataObj$usedObj[[MDS_TYPE]] )
				#browser()
				tab <- t(as.matrix(dataObj$usedObj[[MDS_TYPE]][[clusterby]]))
			}
			else {
				stop( paste("Sorry, the data type onwhat = ",onwhat," is not supported", sep="'") )
			}
			
			if ( ! is.null(useGrouping) ) {
				clusters <- dataObj$samples[,useGrouping]
				if ( is.factor( clusters)) {
					clusters = as.numeric(clusters)
				}
				dataObj <- colors_4 (dataObj, useGrouping )
			}

			## here I need to check that I have no var == 0 samples
			prob <- which(apply(tab, 2, var) == 0)
			if ( length(prob) > 0 ){
				for ( id in prob ){
					tab[,id] = stats::runif(nrow(tab), -19, -16)
				}
			}
			if ( ctype=='hierarchical clust'){
				hc <- stats::hclust(as.dist( 1- stats::cor(tab, method='pearson') ),method = cmethod)
				clusters <- stats::cutree(hc,k=groups.n)
			}else if (  ctype=='kmeans' ) {
				hc <- stats::hclust(as.dist( 1- stats::cor(tab, method='pearson') ),method = cmethod)
				clusters <- stats::kmeans( t(tab) ,centers=groups.n)$cluster
			}else if ( ctype =='mclust' ) {
				
				hc <- mclust::hc( stats::as.dist( 1- stats::cor(tab, method='pearson') ) )
				clusters <- mclust::hclass(hc, groups.n)
			}
			else { stop( paste('ctype',ctype, 'unknown!' ) )}
			
			if ( is.null(useGrouping) ){
				## define the group name n and populate the samples table
				if ( is.null(name)){
					if(is.null(dataObj$usedObj[['auto_clusters']])){
						dataObj$usedObj[['auto_clusters']] = 0
					}
					dataObj$usedObj[['auto_clusters']] <- dataObj$usedObj[['auto_clusters']] +1
					name <- paste( 'auto_clusters', 
							dataObj$usedObj[['auto_clusters']] ,sep='.')
				}
				dataObj$samples <- cbind ( dataObj$samples, clusters )
				colnames(dataObj$samples)[ncol(dataObj$samples)] = name
				clusters <- dataObj$usedObj[['clusters']]
				dataObj$usedObj$usedGrouping <- name
				dataObj <- colors_4(dataObj, name )
				print ("used a new grouing")
			}else {
				print ( "reusing old grouping" )
				dataObj$usedObj$usedGrouping <- useGrouping
			}
			## now I want to create some gene clusters too based on hclust only
#			if ( is.null(dataObj$annotation$'hclust Order')){
#				hcG <- hclust(as.dist( 1- cor(dataObj$data(), method='pearson') ),method = cmethod )
#				dataObj$annotation$'hclust Order' <- hcG$order
#				dataObj$annotation$'hclust 5 groups' <- factor(cutree(hcG,k=5) )
#				dataObj$annotation$'hclust 10 groups' <- factor(cutree(hcG,k=10) )
#				for ( i in c('hclust Order', 'hclust 5 groups', 'hclust 10 groups' )){
#					dataObj <- colors_4(dataObj, i )
#				}
#			}
			dataObj$usedObj[['clusters']] <- clusters
			dataObj$usedObj[['hc']] <- hc
			invisible(dataObj)
} )
