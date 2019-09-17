#' @name createDuplicates
#' @aliases createDuplicates,BioData-method
#' @rdname createDuplicates-methods
#' @docType methods
#' @description Create artificial duplicates using the raw data stored in the object.
#' @param x a SingleCells object
#' @param group which group should be the basis to create duplicates (only inter group duplicates will be created)
#' @param randomCellsFrac how many random cell should be created default=.1
#' @title description of function createDuplicates
#' @export 
if ( ! isGeneric('createDuplicates') ){setGeneric('createDuplicates', ## Name
	function ( x, group, randomCellsFrac=.1 ) { 
		standardGeneric('createDuplicates')
	}
) }

setMethod('createDuplicates', signature = c ('BioData'),
	definition = function ( x, group, randomCellsFrac=.1 ) {
	if ( !is.null(x$raw) ){
		tab = x$raw
	}else {
		tab = x$dat
	}
	## now get the raw data and merge
	total = table( x$samples[,group])
	total = round(total * randomCellsFrac)
	res = NULL
	groupCombo = NULL
	for ( n in levels( x$samples[,group])) {
		THIS= match( x$samples[, group], n )
		OTHER = which(is.na(THIS))
		THIS=which(is.na(THIS) == F)
		sample1 = sample( THIS, total[n] )
		sample2 = sample( OTHER, total[n] )
		add = tab[,sample1] +  tab[,sample2]
		colnames(add) = paste(sep='_', colnames(tab)[sample1], colnames(tab[sample2]))
		groupCombo= c( groupCombo, paste(sep=" + ", as.vector(x$samples[sample1,group]), as.vector(x$samples[sample2,group])  ))
		res= Matrix::cBind( res, add)
	}
	rownames(res) = rownames( tab )
	merged <- as_BioData(matrix( 1, ncol=2, nrow=2, dimnames= list( c('A','B'), c('C','D')) ))
	merged$dat = res
	merged$annotation= data.frame( 'GeneID' = rownames(res), 'info'=rep( 'ArtificialDuplicate', nrow(res)) )
	
	merged$samples=data.frame( 'sampleName' = colnames(res), "Group_combo" = groupCombo)
	class(merged) = class(x)
	merged
} )
