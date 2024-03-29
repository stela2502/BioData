require('R6')

#' Class a simple interface to biological data (numeric) and rich annotation for both columns (samples) and rows (values)
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom grDevices dev.off png rainbow x11
#' @importFrom graphics legend axis barplot image par pie abline hist layout lines mtext plot plot.new rect text title
#' @importFrom stats quantile as.dendrogram density dist hclust median order.dendrogram reorder sd
#' @export
#' @keywords BioData
#' @return Object of \code{\link{R6Class}} to store BioData.
#' @format \code{\link{R6Class}} object.
#' @examples
#' set.seed(1)
#' dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
#' colnames(dat) <- paste('Sample', 1:10)
#' rownames(dat) <- paste( 'gene', 1:100)
#' samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
#' annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )
#' x <- BioData$new( cbind(annotation,dat), 
#' 	Samples=samples, name="testObject",namecol='sname', outpath = "" )
#' @field data the numerical data as sparse Matrix
#' @field raw the raw numerical data as sparse Matrix
#' @field raw_intronic the raw numerical intronic data as sparse Matrix
#' @filed zscored the zscored data as sparse matrix (mean 10.0 sd 1.0)
#' @field samples the sample annotation as data.frame
#' @field annotation the row annotation as data.frame
#' @field usedObj a multi purpose list to store whichever ananlyis results do not fit in the stats list
#' @field stats all stats with one result for each data row
#' @field ranks currently unused
#' @field logged a logical value stating if the data has been logged
#' @field snorm a logical value stating if the object has been sample or cell normalized
#' @field rownamescol which column in the annotation table contains the data rownames
#' @field sampleNamesCol which column in the samples table contains the data rownames
#' @field outpath where to store the object
#' @field name the name of the object (e.g. project name)
#' @field drop=c('MDS') currently unused
#' @field version the version string of the BioData library that created this object
#' @export 
BioData <- #withFormalClass(
		R6::R6Class(
				'BioData',
				class = TRUE,
				public = list ( 
						#' @field dat the numerical data as sparse Matrix
						dat=NULL,
						#' @field raw the raw numerical data as sparse Matrix
						raw=NULL,
						#' @field raw_intronic the raw numerical intronic data as sparse Matrix
						raw_intronic=NULL,
						#' @filed zscored the zscored data as sparse matrix (mean 10.0 sd 1.0)
						zscored=NULL,
						#' @field samples the sample annotation as data.frame
						samples=NULL,
						#' @field annotation the row annotation as data.frame
						annotation=NULL,
						#' @field ranks currently unused
						ranks=NULL,
						#' @field logged a logical value stating if the data has been logged
						logged=FALSE,
						#' @field stats all stats with one result for each data row
						stats=NULL,
						#' @field snorm a logical value stating if the object has been sample or cell normalized
						snorm=FALSE,
						#' @field snorm a logical value stating if the object has been sample or cell normalized
						rownamescol=NULL,
						#' @field rownamescol which column in the annotation table contains the data rownames
						sampleNamesCol=NULL,
						#' @field outpath where to store the object
						outpath='../outpath/',
						#' @field name the name of the object (e.g. project name)
						name='BioData',
						#' @field usedObj a multi purpose list that storeas everything we forgot a slot for
						usedObj = NULL,
						drop=c('MDS'),
						version=NULL,
						#' @description
						#' print a summary of the object
						print = function () {
							cat (paste("An object of class", paste(collapse="::",rev(class(self))),"\n" ) )
							cat("named ",self$name,"\n")
							cat (paste( 'with',nrow(self$dat),'genes and', ncol(self$dat),' samples.'),"\n")
							if ( ! is.null(self$raw) ){
								cat ( "raw ")
								if ( !is.null(self$zscored)) {
									cat( " - and z.scored")
								}
								cat (" data is also stored\n")
							}
							else if ( ! is.null(self$zscored)){
								cat( "z.scored data is also stored\n")
							}
							cat (paste("Annotation datasets (",paste(dim(self$annotation),collapse=','),"): '",paste( colnames(self$annotation ), collapse="', '"),"'  ",sep='' ),"\n")
							cat (paste("Sample annotation (",paste(dim(self$samples),collapse=','),"): '",paste( colnames(self$samples ), collapse="', '"),"'  ",sep='' ),"\n")
							if ( length(names(self$stats)) > 0 ){
								cat ( "P values were calculated for ", length(names(self$stats)) -1, " condition(s)\n")
							}
							MDSnames = NULL
							for( listID in grep('^MDS', names(self$usedObj)) ){
								if ( length(  names(self$usedObj[[listID]])) > 0) {
									cat (paste( sep="",
										names(self$usedObj)[listID], 
										" entries: '",
										paste( MDSnames,collapse="', '"), 
										"'") 
									)
								}
							}
						},
						#' @description
						#' initialize the object - depricated - use the as_BioData functions instead
						#' @param dat the expression sparse Matrix
						#' @param Samples a data.frame describing the column data
						#' @param annotation a data.frame describing the row data
						#' @param name the name of the project - make it phony - will be used as filename
						#' @param namecol the column in the Samples data that is the rownames of the data Matrix
						#' @param namerow the column in the annotation data that is the rownames of the data Matrix
						#' @param outpath the path the file is stored in (defaults to getwd())

						initialize = function (dat,Samples, annotation=NULL, name='BioData', namecol=NULL, namerow= 'GeneID', outpath = ''  ){
							
							S <- Samples
							
							if ( is.null(namecol)){
								stop("Please specify the name of the sample name column (namecol)")
							}
							n <- make.names(as.vector(S[,namecol]))
							if ( is.null(annotation)) {
								mat <- match( as.vector(S[,namecol]), colnames(dat))
								if ( sum(is.na(mat)) > 0 ) {
									stop(paste( 'The samples', 
													paste( as.vector(S[,namecol])[is.na(mat)], collapse=', '),
													'Do not have a data column in the "dat" data.frame' ) 
									)
								}
							
								self$dat =  dat[, mat ]
								annotation <- dat[, is.na(match( colnames(dat), as.vector(S[,namecol]) ))==T ]
							}else {
								self$dat =  dat
							}
							
							if ( class(annotation) == 'factor'){
								annotation <- data.frame( annotation )
								colnames(annotation) <- namerow
							}
							if ( class(annotation) == 'character'){
								annotation <- data.frame( annotation )
								colnames(annotation) <- namerow
							}
							
							if ( outpath == '' ){
								outpath = getwd()
							}
							if ( ! file.exists(outpath)){
								dir.create( outpath )
							}
							
							if ( is.null(dim(annotation))){
								## this xcould be a problem... hope we at least have a vector
								if ( length(annotation) == nrow(self$dat)) {
									rownames(self$dat) <- annotation
								}
								else {
									stop ( "Sorry, please cbind the rownames to the expression values before creating this object!")
								}
							}else{
								rownames(self$dat) <- annotation[,namerow]
							}
							self$samples <- S
							self$name <- name
							self$annotation <- annotation
							self$rownamescol <- namerow
							self$outpath <- outpath
							self$usedObj <- list()
							colnames(self$dat) <- make.names(self$forceAbsoluteUniqueSample ( as.vector(S[, namecol]) ))
							self$samples$samples <- colnames(self$dat)
							
							if ( class(self$dat) != 'dgCMatrix' ) {
								self$dat <- Matrix::Matrix( as.matrix(self$dat), sparse=T ) ## should save up to 80% of memory!
								rm(dat)
								gc()
							}
							
							self$sampleNamesCol <- 'samples'
							self$version = utils::sessionInfo('BioData')$otherPkgs$BioData$Version
							self$force.numeric()
						},
						#' reorder a mds object based on an ordering ID
						#' This function should never be called by the user!!!
						#' @docType methods
						#' @param mdsName the name of the MDS
						#' @param ids the new order
						#' @param mdsClass the MDS class name like MDS_PCA100
						reorder.mds = function( mdsName, ids, mdsClass ){
							self$usedObj[[mdsClass]][[mdsName]] = self$usedObj[[mdsClass]][[mdsName]][ids,]
						},
						#' reorder a BioData object based on annotation information.
						#' If the ordering data is all numeric the numeric values and not the factor levels will be used!
						#' @aliases reorder.genes,BioData-method
						#' @rdname reorder.genes-methods
						#' @docType methods
						#' @description this function reorderes the BioData object based on a column in the annotation table (e.g. for plotting)
						#' @param column the annotation column to reorder on
						reorder.genes = function (column) {
							suppressWarnings( {
							if (is.factor(self$annotation[,column]) &
								all( is.null(as.numeric(as.vector(self$annotation[,column])))) == FALSE ){
								idx <-  order(as.numeric(as.vector(self$annotation[,column])))
							}else {
								idx <-  order( self$annotation[,column])
							}
							})
							self$dat <- self$dat[ idx ,]
							if( !is.null(self$zscored)) {
								self$zscored <- self$zscored[ idx, ]
							}
							if( !is.null(self$raw)) {
								self$raw <- self$raw[ idx, ]
							}
							t <- self$annotation[ idx,]
							if ( class(t) == 'factor' ) {
								t <- data.frame(t)
								colnames(t) <- colnames(self$annotation)
								}
							if ( ! is.null(self$stats)) {
								lapply( names(self$stats), function(n) { self$stats[[n]] <-self$stats[[n]][idx,] })
							}
							for ( n in names(self$usedObj$MDSgene)){
								self$reorder.mds(n, idx,'MDSgene')
							}
							for ( n in names(self$usedObj$MDSgenes_PCA100)){
								self$reorder.mds(n, idx,'MDSgenes_PCA100')
							}
							if ( !is.null(self$usedObj$prGenes) ) {
								self$usedObj$pr = self$usedObj$prGenes[idx,]
							}
							self$annotation <- t
						},
						#' reorder a BioData object based on sample information.
						#' If the ordering data is all numeric the numeric values and not the factor levels will be used!
						#' @aliases reorder.samples,BioData-method
						#' @rdname reorder.samples-methods
						#' @docType methods
						#' @description this function reorderes the BioData object based on a column in the samples table (e.g. for plotting)
						#' @param column the samples column to reorder on
						reorder.samples = function(column) {
							suppressWarnings( {
							if (is.factor(self$samples[,column]) &
								all( is.null(as.numeric(as.vector(self$samples[,column])))) == FALSE ){
								idx <-  order(as.numeric(as.vector(self$samples[,column])))
							}else {
								idx <-  order( self$samples[,column])
							}
							})
							self$dat <- self$dat[ , idx ]
							if( !is.null(self$zscored)) {
								self$zscored <- self$zscored[ , idx]
							}
							if( !is.null(self$raw)) {
								self$raw <- self$raw[ , idx]
							}
							for ( n in names(self$usedObj$MDS)){
								self$reorder.mds(n, idx,'MDS')
							}
							for ( n in names(self$usedObj$MDS_PCA100)){
								self$reorder.mds(n, idx,'MDS_PCA100')
							}
							if (isS4(self$usedObj$pr) ) {
								self$usedObj$pr@scores = self$usedObj$pr@scores[idx,]
							}
							else if ( !is.null(self$usedObj$pr$x) ) {
								self$usedObj$pr$x = self$usedObj$pr$x[idx,]
							}
							self$samples <- self$samples[idx ,]
						},
						#' @description
						#' data accessor function either reports the dat spot or the zscored if existing
						data = function(){
							if ( is.null(self$zscored)){
								self$dat
							}else {
								self$zscored
							}
						},
						#' @description
						#' data accessor function either reports the dat spot or the raw if existing
						rawData = function(){
							if ( is.null(self$raw) ) {
								self$dat
							} else {
								self$raw
							}
						},
						#' @description
						#' useless and depricated
						force.numeric = function() {
							## useless for a Matrix - that one can only be numeric ;-)
							#lapply(colnames(self$dat), function(x) 
							#		{
							#			self$dat[,x] = as.numeric(as.character(self$dat[,x]))
							#		})
							self
						},
						#' @description useless
						pwd = function () {
							system( 'pwd > __pwd' )
							t <- utils::read.delim( file = '__pwd', header=F)
							t <- as.vector(t[1,1])
							t <- paste(t,"/",sep='')
							unlink( '__pwd')
							t
						},
						#' @description 
						#' likely useless and sotherwise implemented, but makes sure all entries in the vector are unique.
						#' by default adds _1, _2, ... , _n to the not unique names
						#' @param x the vector
						#' @param separator by default '_'
						forceAbsoluteUniqueSample = function( x ,separator='_') {
							last = ''
							ret <- vector(length=length(x))
							for ( i in 1:length(x) ){
								if ( is.null(ret) ){
									last = x[i]
									ret[i] <- last
								}
								else{
									last = x[i]
									if ( ! is.na(match( last, ret )) ){
										last <- paste(last,separator,sum( ! is.na(match( x[1:i], last )))-1, sep = '')
									}
									ret[i] <- last
								}
							}
							ret
						}
				)
		)

#' Class interface for the MINT output
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords MINT tRNA
#' @return Object of \code{\link{R6Class}} to store MINT results.
#' @format \code{\link{R6Class}} object.
#' @examples
#' set.seed(1)
#' dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
#' colnames(dat) <- paste('Sample', 1:10)
#' rownames(dat) <- paste( 'gene', 1:100)
#' samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
#' annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )
#' x <- tRNAMINT$new( cbind(annotation,dat), 
#'   Samples=samples, name="testObject",namecol='sname', outpath = "" )
#' @field data the numerical data as data.frame
#' @field samples the sample annotation as data.frame
#' @field annotation the row annotation as data.frame
#' @field usedObj a multi purpose list to store whichever ananlyis results do not fit in the stats list
#' @field stats all stats with one result for each data row
#' @export tRNAMINT
tRNAMINT <-
		R6::R6Class( 'tRNAMINT',
				inherit = BioData,
				class = TRUE
		)

#' Class interface for single cell data
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords single cell, NGS
#' @return Object of \code{\link{R6Class}} to store single cell data.
#' @format \code{\link{R6Class}} object.
#' @examples
#' set.seed(1)
#' dat = matrix(rnorm(1000),ncol=10)
#' dat = round(dat)
#' dat = dat - min(dat)
#' dat = data.frame( dat ) 
#' colnames(dat) <- paste('Sample', 1:10)
#' rownames(dat) <- paste( 'gene', 1:100)
#' samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
#' annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )
#' x <- SingleCells$new( cbind(annotation,dat), 
#'   Samples=samples, name="testObject",namecol='sname', outpath = "" )
#' @field data the numerical data as data.frame
#' @field samples the sample annotation as data.frame
#' @field annotation the row annotation as data.frame
#' @field usedObj a multi purpose list to store whichever ananlyis results do not fit in the stats list
#' @field stats all stats with one result for each data row
#' @export SingleCells
SingleCells <-
		R6::R6Class( 'SingleCells',
				inherit = BioData,
				class = TRUE,
				public = list ( 
						force.numeric = function(...) {
							## check if all data is int
							if ( length(which(apply(self$dat,1,function(a) { all.equal(a, round(a)) == T } ) == F )) != 0 ){
								stop( "A SingleCells object can only be created from raw count data (integers).")
							}
							super$force.numeric()
						}
				)
		)


#' Class interface for two coor microarray data
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords Affy MicroArray
#' @return Object of \code{\link{R6Class}} to store two color microarray data (gene level)
#' @format \code{\link{R6Class}} object.
#' @examples
#' set.seed(1)
#' dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
#' colnames(dat) <- paste('Sample', 1:10)
#' rownames(dat) <- paste( 'gene', 1:100)
#' samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
#' annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )
#' x <- MicroArray$new( cbind(annotation,dat), 
#'   Samples=samples, name="testObject",namecol='sname', outpath = "" )
#' @field data the numerical data as data.frame
#' @field samples the sample annotation as data.frame
#' @field annotation the row annotation as data.frame
#' @field usedObj a multi purpose list to store whichever ananlyis results do not fit in the stats list
#' @field stats all stats with one result for each data row
#' @export MicroArray
MicroArray <-
		R6::R6Class( 'MicroArray',
				inherit = BioData,
				class = TRUE
		)

## obtained from https://rappster.wordpress.com/2015/04/03/r6s3-and-s4-getting-the-best-of-both-worlds/

.onAttach <- function(libname, pkgname) {
	packageStartupMessage(paste("Loading BioData version", packageVersion('BioData') ))
	where <- as.environment("package:BioData")
	clss <- list(
			c("BioData", "R6"),
			c("SingleCells", "BioData"),
			c("tRNAMINT", "BioData"),
			c("MicroArray", "BioData")
	)
	## Ensure clean initial state for subsequent package loads
	## while developing //
	sapply(clss, function(cls) {
				idx <- sapply(cls, isClass )
				suppressWarnings(try(sapply(cls[idx], removeClass,
										where = where), silent = TRUE))
			})
	## Set formal class equivalent //
	sapply(clss, function(cls) {
				try(setOldClass(cls, where = where), silent = TRUE)
			})
#	r6x::formalizeClasses()
}


#try(t <- BioData$new(),silent=T)