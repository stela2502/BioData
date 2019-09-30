#' Based on a PCA representation of the BioData object this function runs 'iter'
#' kmeans clusterings and collects the mean fraction of 'cells' overlapping with each other cell.
#' 
#' The function is returung this fraction of overlap for each cells and the higher the values the better the overlap.
#' 
#' @name findOverlapWithCells
#' @aliases findOverlapWithCells,BioData-method
#' @rdname findOverlapWithCells-methods
#' @docType methods
#' @description PCA and kmeans based similarity finder  
#' @param x the BioData obejct
#' @param cells the cells to find similarity for (do not take too view)
#' @param cname the column name to store the similarity values in default="SimilarityToTargetCells"
#' @param iter the amount of kmean runs for this analysis default=100
#' @param group.n the amount of kmeans groups to be created
#' @title description of function findOverlapWithCells
#' @export 
if ( ! isGeneric('findOverlapWithCells') ){setGeneric('findOverlapWithCells', ## Name
	function ( x, cells, cname="SimilarityToTargetCells", iter=100, groups.n=20 ) { 
		standardGeneric('findOverlapWithCells')
	}
) }

setMethod('findOverlapWithCells', signature = c ('BioData'),
	definition = function ( x, cells, cname="SimilarityToTargetCells", iter=100, groups.n=20  ) {
		
		pb <- progress_bar$new(total = iter)
		
		ids = match(cells, colnames(x))
		if ( length(which(is.na(ids))) > 0 ){
			stop( paste( "The cells", paste( sep=", ", cells[which(is.na(ids))]), "are missing in the object"))
		}
		old_ids = grep("PCA kmeans similarity", colnames(x$samples) )
		if ( length( old_ids) > 0 ) {
			x$samples = x$samples[,- old_ids]
		}
		res = matrix(0, ncol=iter, nrow=ncol(x$dat))
		for ( i in 1:100 ) {
			pb$tick()
			Cname = paste("PCA kmeans similarity",i)
			suppressMessages ( { ## repress print out
			clusters( x , clusterby = "PCA", groups.n = 20, ctype = "kmeans", onwhat= 'MDS', name = Cname )
				}
			)
			cluster = x$samples[,Cname]
			t1 = table(cluster[-ids])
			t2 = table(cluster[ids])
			if (! all(table(c(names(t2), names(t1)) )==2) ) {
				tmp = t1
				m = match( names(t1), names(t2))
				for ( a in 1:length(t1)) {
					if ( is.na(m[a]) ){
						tmp[a] = 0
					}else {
						tmp[a]=t2[m[a]]
					}
				}
				t2 = tmp
			}
			group_stat = t2/length(ids)
#       browser()
			m = match(cluster, names(group_stat))
			res[, i] = res[, i] + group_stat[m]
		}
		stat = apply( res, 1, function(d) { sum(d) / 100 } )
		stat[which(is.na(stat))] = 0
		x$samples[,cname] = stat
		stat
		
} )
