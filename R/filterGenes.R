#' Mitocondrial genes are expressed not from the diploid genome, 
#' but from multiple copies of the mitochondrial genome.
#' A similar problem arises with the ribosomal genes which are 
#' expressed from multiple locis in the diploig genome and therfore also spike at much higher 
#' levels than the diploid genes.
#' 
#' I am not 100% sure about the ribosomal genes, but the mitochondial ones should definitely be filtered out!
#' 
#' @name filterGenes
#' @aliases filterGenes,BioData-method
#' @rdname filterGenes-methods
#' @docType methods
#' @description Filter out genes following a different expression profile
#' @param x the BioData object
#' @param genes a list of pattern matches to the egnes you want to remove
#' @title description of function filterGenes
#' @examples 
#' \dontrun{
#' ## merged is a BioData object with mouse data
#' filterGenes( merged, c('^mt-', 'Rp[ls]'))
#' ## merged is a BioData object with human data
#' filterGenes( merged, c('^MT', 'RP[LS]'))
#' }
#' @export 
if ( ! isGeneric('filterGenes') ){setGeneric('filterGenes', ## Name
	function (x, genes ) { 
		standardGeneric('filterGenes')
	}
) }

setMethod('filterGenes', signature = c ('BioData'),
	definition = function (x, genes ) {
	
	x$samples$Origcounts <- FastWilcoxTest::ColNotZero( x$dat)
	filtered= 0
	for ( gene in genes ) {
		specCount = paste(gene,'counts', sep="_")
		specFraction = paste(gene,'fraction', sep="_")
		mtgenes = grep( gene, rownames(x))
		filteres = filtered + length(mtgenes)
		x$samples[, specCount] <- FastWilcoxTest::ColNotZero( x$dat[mtgenes,])
		x$samples[, specFraction] = x$samples[, specCount] / x$samples$Origcounts
	}
	for ( gene in genes ){
		reduceTo(x, what='row', to= rownames(x$dat) [ - grep (gene, rownames(x$dat)) ] , name=paste( x$name, "Rem", gene,sep="_") )
	}
	x$samples$counts <- FastWilcoxTest::ColNotZero( x$dat)
	print ( paste( filtered ,"genes removed from the dataset"))
	invisible(x)
} )
