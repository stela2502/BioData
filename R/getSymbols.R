#' @name getSymbols
#' @aliases getSymbols,BioData-method
#' @rdname getSymbols-methods
#' @docType methods
#' @description Function copied from https://github.com/PMBio/scLVM/blob/master/R/scLVM/R/util.R
#' is a simple helper to get Gene Symbols from ENSEMBL ids using either org.Mm.egENSEMBL2EG or org.Hs.egENSEMBL2EG, both BioConductor packages.
#' @param ensIds the ids to check
#' @param species the organism to use either default= 'mMus' or 'Hs'
#' @title description of function getSymbols
#' @export 
setGeneric('getSymbols', ## Name
	function (ensIds, species = 'mMus') { ## Argumente der generischen Funktion
		standardGeneric('getSymbols') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('getSymbols', signature = c ('character'),
	definition = function (ensIds, species = 'mMus') {
	if(!(species %in%c('mMus','Hs'))){stop("'species' needs to be either 'mMus' or 'Hs'")}
	
	if(species=='mMus'){
		require(org.Mm.eg.db)
		x <- org.Mm.egSYMBOL
		xxenseg <- AnnotationDbi::as.list(org.Mm.egENSEMBL2EG)}else{
		require(org.Hs.eg.db)
		x <- org.Hs.egSYMBOL
		xxenseg <- AnnotationDbi::as.list(org.Hs.egENSEMBL2EG)        
	}
	# Get the gene symbol that are mapped to an entrez gene identifiers
	gene_names = ensIds
	
	mapped_genes <- mappedkeys(x)
	# Convert to a list
	xx <- as.list(x[mapped_genes])  
	gene_syms=unlist(xx[unlist(xxenseg[gene_names])])
	gene_names_list<-(lapply(xxenseg[gene_names],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
	sym_names=unlist(lapply(xx[unlist(gene_names_list)],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
	sym_names[is.na(sym_names)]=gene_names[is.na(sym_names)]
	
	as.vector(sym_names)
} )
