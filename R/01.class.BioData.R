require('R6')
#require("r6x")

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
#' x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )
#' @field data the numerical data as data.frame
#' @field samples the sample annotation as data.frame
#' @field annotation the row annotation as data.frame
#' @field usedObj a multi purpose list to store whichever ananlyis results do not fit in the stats list
#' @field stats all stats with one result for each data row
#' @export 
BioData <- #withFormalClass(
		R6Class(
		'BioData',
		class = TRUE,
		public = list ( 
				data=NULL,
				raw=NULL,
				samples=NULL,
				annotation=NULL,
				ranks=NULL,
				logged=FALSE,
				stats=NULL,
				snorm=FALSE,
				zscored=FALSE,
				rownamescol=NULL,
				sampleNamesCol=NULL,
				outpath='../outpath/',
				name='BioData',
				usedObj = NULL,
				drop=c('MDS'),
				print = function (...) {
					cat (paste("An object of class", paste(collapse="::",rev(class(self))),"\n" ) )
					cat("named ",self$name,"\n")
					cat (paste( 'with',nrow(self$data),'genes and', ncol(self$data),' samples.'),"\n")
					cat (paste("Annotation datasets (",paste(dim(self$annotation),collapse=','),"): '",paste( colnames(self$annotation ), collapse="', '"),"'  ",sep='' ),"\n")
					cat (paste("Sample annotation (",paste(dim(self$samples),collapse=','),"): '",paste( colnames(self$samples ), collapse="', '"),"'  ",sep='' ),"\n")
					if ( length(names(self$stats)) > 0 ){
						cat ( "P values were calculated for ", length(names(self$stats)) -1, " condition(s)\n")
					}
				},
				initialize = function (dat,Samples, name='BioData', namecol=NULL, namerow= 'GeneID', outpath = ''  ){
					
					S <- Samples
					
					if ( is.null(namecol)){
						stop("Please specify the name of the sample name column (namecol)")
					}
					n <- make.names(as.vector(S[,namecol]))
					mat <- match( as.vector(S[,namecol]), colnames(dat))
					if ( sum(is.na(mat)) > 0 ) {
						stop(paste( 'The samples', 
										paste( as.vector(S[,namecol])[is.na(mat)], collapse=', '),
										'Do not have a data column in the "dat" data.frame' ) 
						)
					}
					
					self$data =  dat[, mat ]
					annotation <- dat[, is.na(match( colnames(dat), as.vector(S[,namecol]) ))==T ]

					if ( class(annotation) == 'factor'){
						annotation <- data.frame( annotation )
						colnames(annotation) <- namerow
					}
					if ( class(annotation) == 'character'){
						annotation <- data.frame( annotation )
						colnames(annotation) <- namerow
					}
					
					if ( outpath == '' ){
						outpath = self$pwd()
					}
					if ( ! file.exists(outpath)){
						dir.create( outpath )
					}
					
					if ( is.null(dim(annotation))){
						## this xcould be a problem... hope we at least have a vector
						if ( length(annotation) == nrow(self$data)) {
							rownames(self$data) <- annotation
						}
						else {
							stop ( "Sorry, please cbind the rownames to the expression values before creating this object!")
						}
					}else{
						rownames(self$data) <- annotation[,namerow]
					}
					self$samples <- S
					self$name <- name
					self$annotation <- annotation
					self$rownamescol <- namerow
					self$outpath <- outpath
					self$usedObj <- list()
					colnames(self$data) <- make.names(self$forceAbsoluteUniqueSample ( as.vector(S[, namecol]) ))
					self$samples$samples <- colnames(self$data)
					
					self$sampleNamesCol <- namecol
					self$force.numeric
				},
				pwd = function () {
					system( 'pwd > __pwd' )
					t <- read.delim( file = '__pwd', header=F)
					t <- as.vector(t[1,1])
					t <- paste(t,"/",sep='')
					unlink( '__pwd')
					t
				},
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
#' x <- tRNAMINT$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )
#' @field data the numerical data as data.frame
#' @field samples the sample annotation as data.frame
#' @field annotation the row annotation as data.frame
#' @field usedObj a multi purpose list to store whichever ananlyis results do not fit in the stats list
#' @field stats all stats with one result for each data row
#' @export tRNAMINT
tRNAMINT <-
		R6Class( 'tRNAMINT',
		inherit = BioData,
		class = TRUE
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
#' x <- MicroArray$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )
#' @field data the numerical data as data.frame
#' @field samples the sample annotation as data.frame
#' @field annotation the row annotation as data.frame
#' @field usedObj a multi purpose list to store whichever ananlyis results do not fit in the stats list
#' @field stats all stats with one result for each data row
#' @export MicroArray
MicroArray <-
		R6Class( 'MicroArray',
				inherit = BioData,
				class = TRUE
		)

## obtained from https://rappster.wordpress.com/2015/04/03/r6s3-and-s4-getting-the-best-of-both-worlds/

.onAttach <- function(libname, pkgname) {
	#packageStartupMessage("Welcome to my package BioData")
	where <- as.environment("package:BioData")
	clss <- list(
			c("BioData", "R6"),
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