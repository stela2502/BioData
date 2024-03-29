#' @name Get_PartialTest_Obj
#' @aliases Get_PartialTest_Obj,BioData-method
#' @rdname Get_PartialTest_Obj-methods
#' @docType methods
#' @description Process the 2 group comparison data and return a summary object
#' @param obj The BioData object
#' @param groupA the outer grouping column name
#' @param groupB the inner grouping column name
#' @param k amount of gene groups (default 1 per 20 genes)
#' @param pcut the cutoff p value default=1e-5
#' @param logfc.threshold  Cpp test option default= .1
#' @param minPct Cpp test option default= .1
#' @param debug lets you debug the mds part default false
#' @title description of function Get_PartialTest_Obj
#' @export 
setGeneric('Get_PartialTest_Obj', ## Name
	function ( obj, groupA, groupB,  k= NULL, pcut=1e-5, logfc.threshold= .1, minPct= .1 , debug=F) { ## Argumente der generischen Funktion
		standardGeneric('Get_PartialTest_Obj') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('Get_PartialTest_Obj', signature = c ('BioData'),
	definition = function ( obj, groupA, groupB, k= NULL, pcut=1e-5, logfc.threshold= .1, minPct= .1 , debug=F) {
        stats <- PartialTests( obj, groupA=groupA, groupB=groupB,
                pcut=pcut, logfc.threshold=logfc.threshold, minPct=minPct)
        print ( paste( length(unique(unlist(stats))), "genes selected" ))
		if ( length(unique(unlist(stats)))<50 ) {
			stop( "less than 50 significant genes detected - please lower the cutoff!")
		}
		forPlot <- reduceTo( obj, what='row', to=unique(unlist(lapply(stats, as.vector ))), copy=T, name=paste(obj$name, groupA,groupB, logfc.threshold, minPct , sep="_", pcut ) )
        ## now create the gene groups:
        for( i in 1:length(stats) ){
                forPlot$annotation[,names(stats)[i]] = "No"
                forPlot$annotation[match( stats[[i]], rownames(forPlot)) ,names(stats)[i]]  = 'Yes'
                forPlot$annotation[,names(stats)[i]] = factor(forPlot$annotation[,names(stats)[i]], levels=c('No','Yes') )
                forPlot$usedObj$colorRange[[names(stats)[i]]] = c('gray', forPlot$usedObj$colorRange[[groupA]][i])
        }
        ## now we need to get a usable grouping into the genes...
        ## try the tSNE + kmeans approach
		if ( is.null(k) ){
        	k = round(length(unique(unlist(stats)))/ 20)
		}
        if (k < 10){ k = 10}
        else if ( k > 100) { k = 100 }
        transpose(forPlot)
		if ( debug ){
			browser()
		}
		DimReduction (forPlot, method='irlba' )
        mds(forPlot, mds.type='TSNE_R') 
        clusters( forPlot , clusterby = "TSNE_R", groups.n = k, ctype = "kmeans", onwhat= 'MDS', name = paste("kmeansTSNE_R clusters_row", sep="_", k ) )
        transpose(forPlot)
        forPlot
} )
