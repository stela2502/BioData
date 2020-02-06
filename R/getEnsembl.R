#' @name getEnsembl
#' @aliases getEnsembl,BioData-method
#' @rdname getEnsembl-methods
#' @docType methods
#' @description Function copied from https://github.com/PMBio/scLVM/blob/master/R/scLVM/R/util.R
#' is a simple helper to get GO annotations from either org.Mm.egGO2EG or org.Hs.egGO2EG, both BioConductor packages.
#' @param term the GO accession number to use e.g. 'GO:0007049' CellCycle
#' @param species the organism to use either default= 'mMus' or 'Hs'
#' @title description of function getEnsembl
#' @return A vector of ENSEMBL gene ids
#' @export 
setGeneric('getEnsembl', ## Name
	function (term, species = 'mMus') { ## Argumente der generischen Funktion
		standardGeneric('getEnsembl') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('getEnsembl', signature = c ('character'),
	definition = function (term, species = 'mMus') {
	if(!(species %in%c('mMus','Hs'))){stop("'species' needs to be either 'mMus' or 'Hs'")}
	
	if(species=='mMus'){
		if(require(org.Mm.eg.db)){
			xxGO <- AnnotationDbi::as.list(org.Mm.egGO2EG)
			x <- org.Mm.egENSEMBL}else{
			stop("Install org.Mm.eg.db package for retrieving gene lists from GO")
		}
	}else{
		if(require(org.Hs.eg.db)){
			xxGO <- AnnotationDbi::as.list(org.Hs.egGO2EG)
			x <- org.Hs.egENSEMBL}else{
			stop("Install org.Hs.eg.db package for retrieving gene lists from GO")
		}
	}
	cell_cycleEG <-unlist(xxGO[term])
	#get ENSEMBLE ids
	
	mapped_genes <- mappedkeys(x)
	xxE <- as.list(x[mapped_genes])
	ens_ids_cc<-unlist(xxE[cell_cycleEG])  
	
	as.vector(ens_ids_cc)
} )
