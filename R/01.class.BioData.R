require('R6')
require("r6x")

#' @name BioData
#' @title BioData
#' @description  An R6 class to store numeric data with corresponding sample and annotation data.
#' @slot objects a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot zscored genes are normalized?
#' @slot snorm samples normalized?
#' @slot usedObj here a set of used and probably lateron important objects can be saved. Be very carful using any of them!
#' @export BioData
BioData <- withFormalClass(
		R6Class(
		'BioData',
		class = TRUE,
		public = list ( 
				data=NULL,
				raw=NULL,
				samples=NULL,
				annotation=NULL,
				ranks=NULL,
				stats=NULL,
				snorm=FALSE,
				zscored=FALSE,
				rownamescol=NULL,
				sampleNamesCol=NULL,
				outpath='../outpath/',
				name='BioData',
				usedObj = NULL,
				drop=c('MDS'),
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
				force.numeric = function ( ) {
					for ( i in 1: ncol(self$data) ) {
						if ( ! is.numeric(self$data[,i]) ) {
							self$data[,i] <- as.numeric(self$data[,i])
						}
					}
					invisible(x)
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
))

#' @name tRNAMINT
#' @title tRNAMINT
#' @description  An R6 class to visualize Expression data.
#' @slot objects a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot zscored genes are normalized?
#' @slot snorm samples normalized?
#' @slot usedObj here a set of used and probably lateron important objects can be saved. Be very carful using any of them!
#' @export tRNAMINT
tRNAMINT <-withFormalClass(
		R6Class( 'tRNAMINT',
		inherit = BioData,
		class = TRUE
))

## obtained from https://rappster.wordpress.com/2015/04/03/r6s3-and-s4-getting-the-best-of-both-worlds/

.onAttach <- function(libname, pkgname) {
#	packageStartupMessage("Welcome to my package BioData")
#	where <- as.environment("package:BioData")
#	clss <- list(
#			c("BioData", "R6"),
#			c("tRNAMINT", "BioData")
#	)
#	## Ensure clean initial state for subsequent package loads
#	## while developing //
#	sapply(clss, function(cls) {
#				idx <- sapply(cls, isClass)
#				suppressWarnings(try(sapply(cls[idx], removeClass,
#										where = where), silent = TRUE))
#			})
#	## Set formal class equivalent //
#	sapply(clss, function(cls) {
#				try(setOldClass(cls, where = where), silent = TRUE)
#			})
	r6x::formalizeClasses()
}

try(t <- BioData$new(),silent=T)