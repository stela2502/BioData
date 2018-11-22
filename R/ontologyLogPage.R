#' @name ontologyLogPage
#' @aliases ontologyLogPage,BioData-method
#' @rdname ontologyLogPage-methods
#' @docType methods
#' @description creates the GO analysis for a gene list and puts it into the report.
#' @param x the BioData object
#' @param genes a list of gene symbols (IMPORTANT)
#' @param ontology which GO ontology to choose from (default = "BP") 
#' @param ... unused
#' @title description of function ontologyLogPage
#' @export 
setGeneric('ontologyLogPage', ## Name
		function ( x, genes, ontology = 'BP', topNodes = 10, GOfname= "GOgenes.csv", ... ) { 
			standardGeneric('ontologyLogPage')
		}
)

setMethod('ontologyLogPage', signature = c ('BioData'),
		definition = function ( x, genes, ontology = 'BP',  topNodes = 10, GOfname= "GOgenes.csv", ... ) {
			## process the ontology for this gene list and add one ontology report page
			error = ""
			
			## for this to work as expected you need an up to date pandoc:
			## https://pandoc.org/installing.html
			
			if ( is.null( x$usedObj$GO2genes)){
				if(x$usedObj$specie =='mouse'){
					if(require(org.Mm.eg.db)){
						db <- org.Mm.eg.db}else{
						stop("Install org.Mm.eg.db package for retrieving gene lists from GO")
					}
				}else if ( x$usedObj$specie=='human'){
					if(require(org.Hs.eg.db)){
						db <- org.Hs.eg.db}else{
						stop("Install org.Hs.eg.db package for retrieving gene lists from GO")
					}
				}else {
					error= paste( "The usedObj$specie",  x$usedObj$specie,  "is up to now not supported in the GO reports function" )
				}
				x$usedObj$GO2genes = mapIds(db, keys(db,'GO'), 'SYMBOL', 'GO', multiVals = 'list')
			}
			
			all = is.na(match(rownames(x$data()), genes ))
			names(all) = rownames(x$data())
			all = factor(all)
			tryCatch({  library("topGO", quietly = TRUE) } ,  
					error = function(e) {
						stop(paste("topGO needed for this function to work. Please install it.\n", e),
								call. = FALSE)
					})
			
			
			x$usedObj$analysis = new("topGOdata", ontology = ontology, allGenes=all 
					,geneSel =  function(x) {x} ,  annot = topGO::annFUN.GO2genes, GO2genes= x$usedObj$GO2genes)
			
			
			resultFisher <- topGO::runTest(x$usedObj$analysis, algorithm = "classic", statistic = "fisher")
			resultKS <- topGO::runTest(x$usedObj$analysis, algorithm = "classic", statistic = "ks")
			resultKS.elim <- topGO::runTest(x$usedObj$analysis, algorithm = "elim", statistic = "ks")
			
			allRes <- topGO::GenTable(x$usedObj$analysis, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim,
					orderBy = "elimKS", ranksOf = "classicFisher", topNodes = topNodes)
			GOI_2_genes <- matrix( 1, nrow=10, ncol=2)
			colnames(GOI_2_genes) = c("GO ID", "Mapping Gene List")
			for( i in 1:nrow(allRes) ) {
				GOI_2_genes[i,1] = allRes[i,1]
				GOI_2_genes[i,2] = paste( 
						unlist( lapply(	intersect( genes,x$usedObj$GO2genes[[allRes[i,1]]]),
										rmdLink, link="https://www.genecards.org/cgi-bin/carddisp.pl?gene=", FALSE ))
						, collapse=" "
				)
			}
			if ( !dir.exists( file.path( x$outpath, 'GOtables') ) ){
				dir.create( file.path( x$outpath, 'GOtables') )
			}
			write.table(GOI_2_genes, sep='\t', quote=F, row.names=F, file= file.path( x$outpath, 'tables',  GOfname ) )
			browser()
			for ( i in 1:nrow(allRes) ) {
				allRes[i,1] = rmdLink(allRes[i,1],"http://amigo.geneontology.org/amigo/term/" )
			}
			#allRes = allRes[,-c(4,5)] ## significant and expected columns do not contain info
			
			
			#write.table(allRes, sep='\t', quote=F, row.names=F, file= file.path( x$usedObj$sessionPath, 'tables', filename(c( n, "GOanalysis.csv") ) ) )
			
			rmd = paste( c(
							#		paste( "##", "GO analysis for grouping", x$usedObj$lastGroup  ),
							paste( "### Genes"),
							paste( collapse="", unlist( lapply( genes,  rmdLink, link="https://www.genecards.org/cgi-bin/carddisp.pl?gene=" ))),
							"",
							paste( "The R package topGO was used to create this output table:"),
							" ",
							" ",
							knitr::kable(allRes, caption=paste("GO analysis for grouping", x$usedObj$lastGroup )),
							" ",
							knitr::kable(GOI_2_genes, caption=paste("The genes mapping to get GO ids" ))
					) , sep="\n")
			
			
			c( genes= GOI_2_genes, stats=allRes, 'Rmd'= rmd )
			
		} )


#' @name rmdLink
#' @aliases rmdLink,BioData-method
#' @rdname rmdLink-methods
#' @docType methods
#' @description creates a linke in the structure [<name>](<link><name>){target='blank'}
#' @param name the displayed text
#' @param link the link address
#' @param lineEnd add a line end at the end of every entry (default =T) 
#' @title easily create an Rmd link that opens in a new window.
#' @export 
if ( ! isGeneric('rmdLink') ){setGeneric('rmdLink', ## Name
			function ( name, link, lineEnd = T ) { 
				standardGeneric('rmdLink')
			}
	) }

setMethod('rmdLink', signature = c ('character'),
		definition = function ( name, link, lineEnd = T  ) {
			
			if ( lineEnd ){
				if(substr(Sys.getenv("OS"),1,7) == "Windows") {
					# set Windows newline
					newLine <- "\r\n"
				}
				else {
					# set non-Windows newline
					newLine <- "\n"
				}
				paste( sep="", "[", name,"](",link,name,"){target='blank'}",newLine)
			}else {
				paste( sep="", "[", name,"](",link,name,"){target='blank'}")
			}
		} )
