#' @name getGeneInfo
#' @aliases getGeneInfo,BioData-method
#' @rdname getGeneInfo-methods
#' @docType methods
#' @description getGeneInfo is a poverful function, that uses the annotation packages RSQLite tables to select information.
#' @title description of function getEnsembl
#' @param x the query string or a vector of query strings
#' @param species either 'mMus', 'Hs' or a AnnotationDbi object
#' @param from the column that corresponds to the x values
#' @param from_tab if the column is defined in more thanone table the specififc table name
#' @param what which column do you want to get (only one)
#' @param what_tab analogue to the from_tab the table name containing the colmn if more than one table contains the info 
#' @return A vector of ids
#' @usage libraray(BioData)
#' CellCycle_gene_symbols <- getGeneInfo(
#' 			"GO:0007049", species = 'mMus', from= 'go_id', from_tab='go' 
#' )
#' CellCycle_ensembl_ids <- getGeneInfo(
#' 			"GO:0007049", species = "Hs", from = "go_id", from_tab = "go", 
#' 			what='ensembl_id', what_tab="ensembl")
#' Complement_and_coagulation_cascades_gene_symbols <- getGeneInfo(
#' 			'04610', species = 'Hs', from = 'path_id', what= 'symbol' )
#' @export 
setGeneric('getGeneInfo', ## Name
	function (x, species = 'mMus', from, from_tab=NULL, what='symbol', what_tab=NULL) { ## Argumente der generischen Funktion
		standardGeneric('getGeneInfo') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('getGeneInfo', signature = c ('character'),
	definition = function (x, species = 'mMus', from, from_tab=NULL, what='symbol', what_tab=NULL ) {
	if(!(species %in%c('mMus','Hs'))){stop("'species' needs to be either 'mMus' or 'Hs'")}
	
	if(species=='mMus'){
		if(require(org.Mm.eg.db)){
			db <- org.Mm.eg.db}else{
			stop("Install org.Mm.eg.db package for retrieving gene lists from GO")
		}
	}else if ( species=='Hs'){
		if(require(org.Hs.eg.db)){
			db <- org.Hs.eg.db}else{
			stop("Install org.Hs.eg.db package for retrieving gene lists from GO")
		}
	}else {
		if ( attributes(class(species))$package == "AnnotationDbi" ) {
			db = species
		}
	}
	
	con <- get('conn', env= db@.xData )
	tables <-  sort(dbListTables(con))
	available <- lapply (tables, function (table) {
				colnames(dbGetQuery(con, paste("select * from ",table," limit 1")))
			}
	)
	names(available) <- tables
	
	get_table <- function( what, what_tab ) {
		
	providing <- NULL
	for ( table in tables ) {
		if ( ! is.na( match( what, available[[table]])) ) {
			providing <- c(providing, table)
		}
	}
	if ( length(providing) == 1) {
		result_tab = providing
	}else if ( length(providing) == 0 ) {
		print( paste( "the column", what, "is not part of the database"))
		return (available)
	}
	else {
		if ( ! is.null(what_tab) & ! is.na(match(what_tab, providing)) ) {
			result_tab = what_tab
		}else {
		stop( paste(
						"sorry, your column", 
						what,"is part of more than one table; please select one from this list:", 
						paste(collapse=", ", providing)
		))
		}
	}
	
	result_tab
	}
	
	# dbListTables(con) # list of tables
	tab = get_table(from, from_tab)
	if ( class(tab) == "list") {
		print ( paste("Sorry, the column", from, "could not be identified in only one table") )
		return (available)
	}
	select_id <- paste("select _id from ",tab," where ",from, " IN ('", paste(x, collapse="', '"),"')",  sep="")
	
	tab =  get_table(what, what_tab)
	if ( class(tab) == "list") {
		print ( paste("Sorry, the column", what, "could not be identified in only one table") )
		return (available)
	}
	select_result <- paste( "select", what,"from", get_table(what, what_tab), "where _id IN (",select_id,")")
	
	result = as.vector(t(dbGetQuery(con, select_result)))
	if ( length(result) == 0 ) {
		print ( "something did not work as expected - please check manually" )
		browser()
	}
	result
} )
