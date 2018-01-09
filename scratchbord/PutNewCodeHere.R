getEnsembl <- function(term, species = 'mMus'){
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
	
	ens_ids_cc
}


getSymbols <- function(ensIds, species = 'mMus'){
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
	
	sym_names
}