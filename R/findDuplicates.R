#' In contrast to the publiushed function this will use the re-sampling based normalization the BioData obejct provides.
#' In addition this function uses the BioData internal VarGenes function.
#' 
#' The results are totally untested and a lot of additional work is likely necessary.
#' 
#' @name findDuplicates
#' @aliases findDuplicates,BioData-method
#' @rdname findDuplicates-methods
#' @docType methods
#' @description This function is a re-implementation of the DuplicateFinder logics
#' @param x a SingleCells object
#' @param group a grouing focusing on the differences between groups. Too many groups would be bad here!
#' @param randomCellsFrac how many cells should be created as random duplicates default=.01
#' @param iterations each iteration takes significant amount of time default= 10
#' @title description of function findDuplicates
#' @export 
if ( ! isGeneric('findDuplicates') ){setGeneric('findDuplicates', ## Name
	function ( x, group, randomCellsFrac=.01, iterations = 10 ) { 
		standardGeneric('findDuplicates')
	}
) }

setMethod('findDuplicates', signature = c ('SingleCells'),
	definition = function ( x, group, randomCellsFrac=.01, iterations = 10 ) {
	## x is a SingleCells object and it will be processed internally as I whish!
	## first I need to get the data volume down to the variable genes.
	x$samples$InTeRnAlNoTuSe = 'notUse'
	
	x$samples$total = factor(rep('total', ncol(x$dat)))
	genes = unlist( getGenesExpressedHigherThanExpected( x,'total', n= 1300 ))
	
	Intern = reduceTo( x, what='row', to=genes, name='Intern', copy=T )
	
	normTo = sum(expm1(x$dat[ which( x$dat[,1] > 0 ),1]))
	
	z.score( Intern )
	removeNeg( Intern, 'zscored' )
	Intern$zscored@x = Intern$zscored@x - 10 ## get 0 centered!
	
	pca = irlba::prcomp_irlba ( Matrix::t( Intern$zscored ), n=20, center=F)
	
	## now I need to create an intern data set and a random merged one
	res = matrix( 0, nrow=ncol(Intern$dat), ncol=iterations)
	rownames(res) = colnames(Intern)	
	
	for (i in 1:iterations){
		print ( paste("run duplicate search", i))
		duplicates <- createDuplicates( x=x, group=group, randomCellsFrac=randomCellsFrac )
		normalize( duplicates, reads=normTo )
		z.score( duplicates )
		removeNeg( duplicates )
		reduceTo( duplicates, what='row', to=genes, name="duplicates", copy=F)
	
		additions <- Matrix::t( duplicates$zscored) %*% pca$rotation
	
		#rownames( pca$x ) = rownames( additions )
		#rownames( additions ) = paste(rownames( additions ), 'copy')
		for (id in 1:20 ){
			cluster =stats::kmeans( as.matrix( rbind( pca$x, additions ) ) ,centers=100)$cluster
			t1 = table(cluster[1:ncol(Intern$dat)])
			t2 = table(cluster[ncol(Intern$dat)+1: length(cluster)])
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
			group_stat = t2/ncol(duplicates$dat)
			m = match(cluster[1:ncol(Intern$dat)], names(group_stat))
			res[, i] = res[, i] + group_stat[m]
		}
		res[, i] = res[, i] / 20
		gc()
	}
	## get the stat for the res table
	stat = apply( res, 1, function(d) { sum(d) / iterations } )
	x$samples[, paste('DuplicateSearch', group) ] = stat
	message( paste( sep="",'The column "DuplicateSearch ', group, '" contains the probability to be a duplicate'))
	
	invisible(Intern) 
} )
