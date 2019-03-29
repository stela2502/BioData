#' @name IdentifyMarkerGenes
#' @aliases IdentifyMarkerGenes,BioData-method
#' @rdname IdentifyMarkerGenes-methods
#' @docType methods
#' @description This method uses whichever stats function was selected for this class using the createStats function.
#' But it compares one group versus all other groups to find marker genes for this group only.
#' This function will add length(group) new stats tables.
#' @param x the BioData object
#' @param gname the samples column name to group on.
#' @param settings a list of slurm parameters to use to run a script (optional)
#' @param names a vector of group name values to compare against all others (default NULL == use all)
#' @title description of function IdentifyMarkerGenes
#' @export 
setGeneric('IdentifyMarkerGenes', ## Name
	function ( x, gname, settings=list(), names=NULL ) { 
		standardGeneric('IdentifyMarkerGenes')
	}
)


setMethod('IdentifyMarkerGenes', signature = c ('BioData'),
		definition = function ( x, gname, settings=list(), names=NULL ) {
			x$name <- stringr::str_replace_all( x$name, '\\s+', '_')
			OPATH <- file.path( x$outpath,stringr::str_replace( x$name, '\\s', '_'))
			opath = file.path( OPATH,gname,"RFclust.mp" )
			putScript <- function( n, ofile ) {
				Rdata = paste(n,'RData', sep='.')
				fileConn<-file( ofile )
				writeLines ( c( 'library(BioData)', 
								'library(RFclust.SGE)',
								paste('set.lock("',file.path(opath,Rdata),'")',sep=''),
								paste(sep="",'load("',file.path(opath,'IdentifyMarkerGenes_tmp.RData'),'")' ) ,
								'#reads object x',
								paste(sep="",'IdentifyMarkerGenes( data, "',n,'" )'),
								"stat = data$stats[[1]]",
								paste(sep="",'save(stat, file="',file.path(opath,Rdata),'")'),
								paste('release.lock("',file.path(opath,Rdata),'")',sep='')
						), con=fileConn )
				close(fileConn)
				cmd <- paste('R CMD BATCH --no-save --no-restore --no-readline --max-ppsize=500000 --', ofile )
				x$usedObj$IdentifyMarkerGenes[[n]] <- file.path(opath,Rdata)
				cmd
			} 
			if ( ! dir.exists(OPATH)){
				dir.create( OPATH )
			}
			if ( ! dir.exists(file.path(OPATH, gname )) ){
				dir.create(file.path(OPATH, gname ) )
			}
			if ( ! dir.exists(file.path(OPATH, gname, "RFclust.mp")) ){
				dir.create(file.path(OPATH, gname,"RFclust.mp" ) )
			}
			if ( is.null(names)) {
				names = unique(tmp$samples[,gname])
			}
			if ( length(names(settings)) == 0){
				tmp <- x$clone()
				for ( n in names) {
					if ( n == 'rest') {
						next
					}
					tmp$stats <- NULL
					gc(FALSE)
					new_g <-  paste( 'IdentifyMarkerGenes',  gname, n )
					print (paste( "Processing:", new_g ))
					g <- rep('rest', ncol(x$dat) )
					g[which(tmp$samples[,gname] == n )] = n
					tmp$samples[,new_g] <- factor( g, levels= c( n, 'rest') ) 
					createStats( tmp, new_g)
					x$stats[[new_g]] = tmp$stats[[1]]
				}
				rm(tmp)
				gc(FALSE)
				return ( invisible(x) )
			}else if ( is.null(x$usedObj$IdentifyMarkerGenes) ) {
				x$usedObj$IdentifyMarkerGenes <- list()
				tmp <- x$clone()
				tmp$stats <- NULL
				groups <- NULL
				rfObj <- RFclust.SGE::RFclust.SGE ( 
						dat=as.data.frame(matrix(0,ncol=10, nrow=10)), 
						SGE=F, slices=1, email="nothing@nowhere.se", tmp.path=opath, 
						name= 'IdentifyMarkerGenes', settings=settings, slurm=T 
				)
				for ( n in names ) {
					new_g <-  paste( 'IdentifyMarkerGenes',  gname, n )
					new_g <- stringr::str_replace_all( new_g, '\\s+', '_')
					print (paste( "Processing:", new_g ))
					g <- rep('rest', ncol(x$dat) )
					g[which(tmp$samples[,gname] == n )] = n
					tmp$samples[,new_g] <- factor( g, levels= c( n, 'rest') ) 
					groups <- c( groups , new_g)
				}
				tmp$outpath = opath
				tmp$name = "IdentifyMarkerGenes_tmp"
				saveObj(tmp)
				#rfObj@debug=TRUE ## for now
				for ( n in groups ) {
					cmd = putScript( n ,file.path(tmp$outpath, paste(sep=".",n,"R") ) )
					RFclust.SGE::writeSLURMscript(rfObj, n , cmd )
					#x$usedObj$IdentifyMarkerGenes[[n]] <- Rdata
				}
			}else { # ! is.null(x$usedObj$IdentifyMarkerGenes) 
				for ( n in names(x$usedObj$IdentifyMarkerGenes)) {
					if ( RFclust.SGE::locked(x$usedObj$IdentifyMarkerGenes[[n]]) ) {
						stop(paste( "Process for grouping", n ,"not finished!" ))
					}
					load(x$usedObj$IdentifyMarkerGenes[[n]])
					x$usedObj$IdentifyMarkerGenes[[n]] <- NULL
					x$stats[[n]] <- stat
				}
			}
			invisible(x)
		} 
)