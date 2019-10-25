MDS_predict <- function( x, mds.type, new, useRaw=FALSE ) {
	## the new data needs to be included into the old one. That is absolutely not simple and needs some re-calculation
	## at the moment this is only possible for PCA and UMAP data as I have not thought about the others...
	if ( mds.type == 'UMAP' ) {
		if ( useRaw ) {
			stop( "useRaw with UMAP has not been implemented - please do that now as you need it ;-)")
		}
		## I need the x$usedObj$pr
		if ( is.null( x$usedObj$pr) ) {
			message( "assuming you use the default DimReduction to create the initial PCA data")
			PCA_name = DimReduction( x )
		}
		
		## and predict into this PCA the new data,
		## extract the new PCA eigenvectors and add them to the old UMAP
		## which I fist need to re-create from the old PCA data
	}
	
}