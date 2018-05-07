runStats_inThread <- function ( x, condition, files=F, A=NULL, B=NULL, covariates=NULL, form=NULL ) {
	## Quite simple - create a script and run it using R CMD Batch &
	ofile_base <- paste( x$name,'_',condition,'runStats_inThread', sep='') 
	fname <- function( name1, ext) { paste( name1, ext, sep=".") }
	
	if ( file.exists(file.path( x$outpath,fname(ofile_base,"pid" ))) ) {
		stop( "The process is still running" )
	}
	else if ( ! file.exists(file.path( x$outpath, fname(ofile_base,"finished" )))) {
		if ( ! file.exists( file.path( x$outpath, fname( x$name, 'RData' ) ) )){
			saveObj(x)
		}
		## create the funciton call
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
		## create the script
		script = paste( sep="\n", "library(BioData)",
				paste( sep="", 'cat(Sys.getpid(),file="',fname(ofile_base,'pid'), ',sep="\n")' ),
				paste( sep="", "data <- loadObj('",x$name,".RData')"  ),
				fcall,
				'saveObj(data)',
				paste( sep="", 'cat(Sys.getpid(),file="', fname(ofile_base,'finished'),'",sep="\n")' ),
				paste( sep="", "unlink('",fname(ofile_base,'pid'),'")' )
		)
		print ("create and execute script")
		cat(script, file=paste(sep='.', ofile_base, 'sh' ) )
		## run the script
		system( 
			paste( 
				"cd", x$outpath,"&&", 
				"R CMD BATCH --no-save --no-restore --no-readline --max-ppsize=500000 --", 
				fname(ofile_base, 'sh') 
			)
		)
	}
	else { ## the script
		print ( "read the data" ) 
		data <- loadObj(file.path( x$outpath, fname( x$name, 'RData' ) ))
		n <- names(data$stats)[length(names(data$stats))]
		m <- match(rownames(data$dat), rownames(x$dat))
		x$stats[[n]] <- data$stats[[n]][m,]
	}
	
	invisible(x)
}