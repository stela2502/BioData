#' @name clusters_gene
#' @aliases clusters_gene,BioData-method
#' @rdname clusters_gene-methods
#' @docType methods
#' @description Culters the data either based on the raw data values or any MDS data type and adds the grouping into the annotation table.
#' @param dataObj the BioData object
#' @param clusterby cluster on raw data or MDS clustered data default="raw"
#' @param useGrouping do nothing and simply use tis grouping default=NULL
#' @param groups.n how many groups should be detected default= 3
#' @param ctype cluster type - either 'hierarchical clust' or 'kmeans' default = 'hierarchical clust'
#' @param onwhat this option has been kept for the Fluidigm data as there FACS data can also be used default = 'Expression'
#' @param cmethod the method to used with the hclust clustering (default = 'ward.D2')
#' @param name the name for the new grouping (default = 'auto_clusters_gene.1:n')
#' @title description of function clusters_gene
#' @export 
if ( ! isGeneric('clusters_gene') ){ setGeneric('clusters_gene', ## Name
	function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3, 
			ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D2', name=NULL ) {## Argumente der generischen Funktion
		standardGeneric('clusters_gene') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'clusters_gene' already defined - no overloading here!")
}

setMethod('clusters_gene', signature = c ('BioData'),
	definition = function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3,
				ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D2', name=NULL ) {
	
			clusters_gene <- NULL
			hc <- NULL
			if(onwhat=="Expression"){
				tab <- t(dataObj$data)
			}
			else {
				stop( paste("Sorry, the data type",mds.type,"is not supported") )
			}
			if ( ! is.null(useGrouping) ) {
				clusters_gene <- dataObj$annotation[,useGrouping]
				if ( is.factor( clusters_gene)) {
					clusters_gene = as.numeric(clusters_gene)
				}
				dataObj <- colors_4 (dataObj, useGrouping )
			}
			else if ( clusterby=="raw"){#...do mds on tab
				if ( ctype=='hierarchical clust'){
					hc <- hclust(as.dist( 1- cor(tab, method='pearson') ),method = cmethod)
					clusters_gene <- cutree(hc,k=groups.n)
				}else if (  ctype=='kmeans' ) {
					hc <- hclust(as.dist( 1- cor(t(tab), method='pearson') ),method = cmethod)
					clusters_gene <- kmeans( dataObj$usedObj[['mds.proj']] ,centers=groups.n)$cluster
				}else if ( ctype =='mclust' ) {
					hc <- hc( as.dist( 1- cor(t(tab), method='pearson') ) )
					clusters_gene <- hclass(hc, 12)
				}
				else { stop( paste('ctype',ctype, 'unknown!' ) )}
			}else { ## now the clusterby is a MDS algorithm name / MDS dataset name
				if ( is.null( dataObj$usedObj$MDSgenes[[clusterby]] ) ) {
					dataObj <- mds( dataObj, onwhat="Expression", mds.type=clusterby, genes=T)
				}
				if ( ctype=='hierarchical clust'){
					hc <- hclust(dist( dataObj$usedObj$MDSgenes[[clusterby]] ),method = cmethod)
					clusters_gene <- cutree(hc,k=groups.n)
				}else if (  ctype=='kmeans' ) {
					hc <- hclust(dist( dataObj$usedObj$MDSgenes[[clusterby]] ),method = cmethod)
					clusters_gene <- kmeans( dataObj$usedObj$MDSgenes[[clusterby]] ,centers=groups.n)$cluster
				}else if ( ctype =='mclust' ) {
					hc <- hc( dataObj$usedObj$MDSgenes[[clusterby]] )
					clusters_gene <- hclass(hc, groups.n)
				}
				else { stop( paste('ctype',ctype, 'unknown!' ) )}
			}
			if ( is.null(useGrouping) ){
				## define the group name n and populate the annotation table
				if ( is.null(name)){
					if(is.null(dataObj$usedObj[['auto_clusters_gene']])){
					dataObj$usedObj[['auto_clusters_gene']] = 0
				}
				dataObj$usedObj[['auto_clusters_gene']] <- dataObj$usedObj[['auto_clusters_gene']] +1
				name <- paste( 'auto_clusters_gene', 
						dataObj$usedObj[['auto_clusters_gene']] ,sep='.')
				}
				dataObj$annotation <- cbind ( dataObj$annotation, clusters_gene )
				colnames(dataObj$annotation)[ncol(dataObj$annotation)] = name
				clusters_gene <- dataObj$usedObj[['clusters_gene']]
				dataObj$usedObj$usedGroupingGene <- name
				dataObj <- colors_4(dataObj, name )
				print ("used a new grouing")
			}else {
				print ( "reusing old grouping" )
				dataObj$usedObj$usedGroupingGene <- useGrouping
			}
			## now I want to create some gene clusters_gene too based on hclust only
#			if ( is.null(dataObj$annotation$'hclust Order')){
#				hcG <- hclust(as.dist( 1- cor(dataObj$data, method='pearson') ),method = cmethod )
#				dataObj$annotation$'hclust Order' <- hcG$order
#				dataObj$annotation$'hclust 5 groups' <- factor(cutree(hcG,k=5) )
#				dataObj$annotation$'hclust 10 groups' <- factor(cutree(hcG,k=10) )
#				for ( i in c('hclust Order', 'hclust 5 groups', 'hclust 10 groups' )){
#					dataObj <- colors_4(dataObj, i )
#				}
#			}
			dataObj$usedObj[['clusters_gene']] <- clusters_gene
			dataObj$usedObj[['hc_gene']] <- hc
			invisible(dataObj)
} )
