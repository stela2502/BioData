#' @name runStats_inThread
#' @aliases runStats_inThread,BioData-method
#' @rdname runStats_inThread-methods
#' @docType methods
#' @description 
#' @param x the BioData object
#' @param condition the samples column to run stats for
#' @param files whether to print the stats output to file ( default =F - depricated)
#' @param A if not all conditions should be used condition A
#' @param B if not all conditions should be used condition B
#' @param covariates should covariates be used - name them here
#' @param form a specific formualr to use? State it here
#' @title description of function runStats_inThread
#' @export 
setGeneric('runStats_inThread', ## Name
	function ( x, condition, files=F, A=NULL, B=NULL, covariates=NULL, form=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('runStats_inThread') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('runStats_inThread', signature = c ('BioData'),
	definition = function ( x, condition, files=F, A=NULL, B=NULL, covariates=NULL, form=NULL ) {
	## Quite simple - create a script and run it using R CMD Batch &
	ofile_base <- paste( x$name,condition, sep='_')
	if ( !is.null(covariates) ) {
		ofile_base <- paste( ofile_base,covariates, collapse="_", sep='_')
	}
	ofile_base <- paste(ofile_base,'_runStats_inThread', sep='_')
	
	ofile_base <- str_replace_all(ofile_base, "\\s+","_") 
	fname <- function( name1, ext) { paste( name1, ext, sep=".") }
	
	if ( file.exists(file.path( x$outpath,fname(ofile_base,"pid" ))) ) {
		stop( "The process is still running" )
	}
	else if ( ! file.exists(file.path( x$outpath, fname(ofile_base,"finished" )))) {
		## create the funciton call
		if ( ! file.exists( file.path( x$outpath, fname( x$name, 'RData' ) ) )){
			saveObj(x)
		}
		fcall <- paste( sep="", "createStats(data, condition='",condition,"'")
		if ( files ) {
			fcall <- paste( sep="", fcall ,', files=T' )
		}
		if (! is.null(A) ) {
			fcall <- paste( sep="", fcall ,', A="',A, '"' )
		}
		if (! is.null(B) ) {
			fcall <- paste( sep="", fcall ,', B="',B, '"' )
		}
		if (! is.null(covariates) ) {
			fcall <- paste( sep="", fcall ,', covariates="',covariates, '"' )
		}
		if (! is.null(form) ) {
			fcall <- paste( sep="", fcall ,', form="',form, '"' )
		}
		fcall <- paste( fcall ,')')
		fcall <- paste( 'eval(',fcall,')')
		## create the script
		script = paste( sep="\n", "library(BioData)",
				paste( sep="", 'cat(Sys.getpid(),file="',fname(ofile_base,'pid'),'")' ),
				paste( sep="", "data <- loadObj('",file.path( x$outpath, fname( x$name, 'RData' ) ),"')"  ),
				'data$stats <- list()',
				fcall,
				'stat_res <- list( name = data$name, stat = data$stats[[1]], stat_name = names(data$stats)[1] )',
				
				paste( sep="", "save(stat_res, file = '", fname(ofile_base,'RData' ),"' )" ),
				paste( sep="", 'cat(Sys.getpid(),file="', fname(ofile_base,'finished'),'")' ),
				paste( sep="", "unlink('",fname(ofile_base,'pid'),"')", "" )
		)
		print (paste ("create and run script", file.path( x$outpath,fname( ofile_base, 'sh' ) ) ) )
		cat(script, file= file.path( x$outpath, fname( ofile_base, 'sh' )) )
		## run the script
		system( 
			paste( sep='',
				"cd ", x$outpath," && ", 
				"R CMD BATCH --no-save --no-restore --no-readline --max-ppsize=500000 -- '", 
				fname(ofile_base, 'sh') , "' &"
			)
		)
	}
	else { ## the script
		print ( "read the data" ) 
		
		load( fname(ofile_base,'RData' ) )
		if ( stat_res$name == x$name ){
			n <- stat_res$stat_name
			m <- match(rownames(stat_res$stat), rownames(x$dat))
			x$stats[[n]] <- stat_res$stat[m,]
		}else {
			stop("The name of the BioData object must not change between stst creation and reading.")
		}
	}
	
	invisible(x)
} )
