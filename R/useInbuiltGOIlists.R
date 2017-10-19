#' @name useInbuiltGOIlists
#' @aliases useInbuiltGOIlists,BioData-method
#' @rdname useInbuiltGOIlists-methods
#' @docType methods
#' @description  An easy function to register the inbuilt (G)enes (O)f (I)nterest lists 'TFs' and
#' @description  'epigenetic' are supported at the moment
#' @param x A BioData object
#' @param name the name of the inbuilt list to use ( either 'TFs' or 'epigenetic' for now)
#' @param gene_col which column in the annotation data contains the Gene.Symbols? 
#' The inbuilt groups are Gene.Symbol based.
#' @title description of function useInbuiltGOIlists
#' @export useInbuiltGOIlists
setGeneric('useInbuiltGOIlists', ## Name
		function (x, name, ...) { ## Argumente der generischen Funktion
			standardGeneric('useInbuiltGOIlists') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('useInbuiltGOIlists', signature = c ('BioData'),
		definition = function (x, name, gene_col=NULL, ... ) {
			
			if ( ! is.na( match(name, colnames(x$annotation)))) {
				stop( "This GIO list has already been defined" )
			}
			
			if ( is.null(gene_col)) {
				stop("Please give me the name of the Gene.Sybol column in the annotation dataset as the inbuilt lists are gene.Symbol based." )
			}
			
			if ( name == "TFs" ){		
				hum_t <- length(which(is.na(match(x$annotation[,gene_col],human.tfs))==F))
				mouse_t <- length(which(is.na(match( x$annotation[,gene_col], mouse.tfs))==F))
				if (hum_t > mouse_t ){
					x = defineGOIs( x, name, human.tfs, gene_col=gene_col )
				}else if ( mouse_t > hum_t ){
					x = defineGOIs( x, name, mouse.tfs, gene_col=gene_col )
				}else {
					stop( "Sorry, but neither inbuilt dataset (Gene Symbols from mouse and humans) do match to the rownames(@data) - please double ckech that.")
				}
				
			}
			else if ( name == 'epigenetic' ) {
				# register 'epigeneic'
				hum_e <- length(which(is.na(match( x$annotation[,gene_col],Epigenetic$HGNC_symbol))==F))
				mouse_e <- length(which(is.na(match( x$annotation[,gene_col],Epigenetic$MGI_symbol ))==F))
				if ( hum_e > mouse_e){
					x = defineGOIs( x, name, Epigenetic$HGNC_symbol, Epigenetic$Target, gene_col=gene_col )
				}else if ( mouse_e > hum_e ){
					x = defineGOIs( x, name, Epigenetic$MGI_symbol, Epigenetic$Target, gene_col=gene_col)
				}else {
					stop( "Sorry, but neither inbuilt dataset (Gene Symbols from mouse and humans) do match to the rownames(@data) - please double ckech that.")
				}
			}
			else if ( name =="CellCycle" ) {
				hum_e <- length(which(is.na(match( x$annotation[,gene_col],CellCycle$Gene.Symbol))==F))
				mouse_e <- length(which(is.na(match( x$annotation[,gene_col],CellCycle$MouseGene ))==F))
				if ( hum_e > mouse_e){
					x = defineGOIs( x, name, CellCycle$Gene.Symbol, CellCycle$X, gene_col=gene_col )
				}else if ( mouse_e > hum_e ){
					x = defineGOIs( x, name, CellCycle$MouseGene, CellCycle$X, gene_col=gene_col )
				}else {
					stop( "Sorry, but neither inbuilt dataset (Gene Symbols from mouse and humans) do match to the rownames(@data) - please double ckech that.")
				}
			}
			else if ( name =="CellSurface" ) {
				hum_e <- length(which(is.na(match( x$annotation[,gene_col], human.CellSurface))==F))
				mouse_e <- length(which(is.na(match( x$annotation[,gene_col], mouse.CellSurface ))==F))
				if ( hum_e > mouse_e){
					x = defineGOIs( x, name, human.CellSurface, gene_col=gene_col )
				}else if ( mouse_e > hum_e ){
					x = defineGOIs( x, name, mouse.CellSurface, gene_col=gene_col )
				}else {
					stop( "Sorry, but neither inbuilt dataset (Gene Symbols from mouse and humans) do match to the rownames(@data) - please double ckech that.")
				}
			}
			else {
				stop ( paste("Sorry, but the gene list", name, "is not defined" ) )
			}
			x
		} 
)
