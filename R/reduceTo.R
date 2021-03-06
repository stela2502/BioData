#' @name reduceTo
#' @aliases reduceTo,BioData-method
#' @rdname reduceTo-methods
#' @docType methods
#' @description The main reduction function can drop both samples and genes using the colnames / rownames of the data tables
#' @param x the NGScollation object
#' @param what reduce to samples or row ids default='row'
#' @param to select these names default=NULL
#' @title description of function reduceTo
#' @export 
if ( ! isGeneric('reduceTo') ){ methods::setGeneric('reduceTo', ## Name
		function ( x, what='row', to=NULL, ... ) { 
			standardGeneric('reduceTo')
		}
)
}else {
	print ("Onload warn generic function 'reduceTo' already defined - no overloading here!")
}

setMethod('reduceTo', signature = c ('BioData'),
		definition = function ( x, what='row', to=NULL, name='reducedTo', copy=FALSE ) {
			if ( copy ) {
				x <- x$clone()
			}
			if ( ! is.null(to)) {
				if ( what =="row") {
					if ( length(which(is.na(match(to,rownames(x$dat)))==F)) > 0 ) {
						useOnly <- match(to, rownames(x$dat))
						not.matched <- to[is.na(useOnly)]
						if ( length(not.matched) > 0 ){
							print (paste('Problematic genes:', paste(not.matched,sep=', ')))
							to <- to[ ! is.na(useOnly)]
							useOnly <- useOnly[ ! is.na(useOnly) ]
						}
						for (n in x$drop){
							if ( ! is.null(x[[n]]) ) {
								x[[n]] <- NULL
							}
							if ( ! is.null(x$usedObj[[n]]) ) {
								x$usedObj[[n]] <- NULL
							}
						}
						x$dat <- x$dat[useOnly,]
						x$annotation <- x$annotation[useOnly,]
						
						if ( ! is.null( x$raw) ) {
							x$raw <- x$raw[useOnly,]
						}
						if ( ! is.null( x$zscored) ) {
							x$zscored <- x$zscored[useOnly,]
						}
						if ( length(x$stats) > 0 ) {
							lapply( names(x$stats) , function(name) {
										if ( length( grep('Cpp_FindAllMarkers' ,name) )== 1 ){
											statIDs <- match(tolower(x$stats[[name]][,'gene']) , tolower(to) )
											statIDs= which(is.na(statIDs) == F)
											add_to_stat( x, x$stats[[name]][statIDs,] , name )
										}else {
											statIDs <- match(tolower(to) ,tolower(x$stats[[name]][,1]) )
											if ( length(which(is.na(statIDs))) > 0 ){
												statIDs <- match(tolower(to) ,tolower(rownames(x$stats[[name]]) ) )
											}
											add_to_stat( x, x$stats[[name]][statIDs,] , name )
										#x$stats[[name]] <- x$stats[[name]][statIDs,]
										}
									} )
						}
						# to be sure that is really OK
						x$usedObj$pr  = x$usedObj$prGenes = NULL
						## but the MDS objects should probably be reduced...
						for ( listIDs in grep( '^MDS', names(x$usedObj)) ) {
							## no, they actually loose there validity and should if anything be manually copied over.
						}
						x$name = name
					}else {
						print (paste( "None of the probesets matched the probesets in",x$name, "-> keep everything!"))
					}
					
					
				}else if ( what =="col" ) {
					#to <- tolower(make.names(to))
					if ( length(which(is.na(match(to,colnames(x$dat)))==F)) > 0 ) {
						useOnly <- match(to, colnames(x$dat))
						not.matched <- to[is.na(useOnly)]
						if ( length(not.matched) > 0 ){
							print (paste('Problematic samples:', paste(not.matched,sep=', ')))
							to <- to[ ! is.na(useOnly)]
							useOnly <- useOnly[ ! is.na(useOnly) ]
						}
						
						x$dat <- x$dat[,useOnly]
						x$samples <- x$samples[useOnly,]
						
						if ( ! is.null( x$raw) ) {
							x$raw <- x$raw[,useOnly]
						}
						if ( ! is.null( x$zscored) ) {
							x$zscored <- x$zscored[,useOnly]
						}
						if ( length(x$stats) > 0 ) {
							x$stats = list()
						}
						for (n in x$drop){
							if ( ! is.null(x[[n]]) ) {
								x[[n]] <- NULL
							}
							if ( ! is.null(x$usedObj[[n]]) ) {
								x$usedObj[[n]] <- NULL
							}
						}
						# to be sure that is really OK
						x$usedObj$pr  = x$usedObj$prGenes = NULL
						## but the MDS objects should probably be reduced...
						## NO this creates a mess! I now have a copy_mds function that can do that if needed! 
						for ( listID in grep( '^MDS', names(x$usedObj)) ) {
							for (n in names(x$usedObj[[listID]])){
								new_name = paste(x$name, sep="_", n)
								## need to modify the old one and then rename!
								try( {	
									x$usedObj[[listID]][[n]] = x$usedObj[[listID]][[n]][useOnly,]
									names(x$usedObj[[listID]])[match(n, names(x$usedObj[[listID]]))] = new_name
								})
							}
						}
						x$name = name
						
					}else {
						print (paste( "None of the names (to) matched the sample names in",x$name, "-> keep everything!"))
						print (paste( "to =",paste(tolower(to), collapse=", ")))
						print (paste( "match to =",paste(tolower(colnames(x$dat)), collapse=", ")))
					}
				}else {
					stop(paste( "the option what='",what,"' is not supported!", sep="") )
				}
			}
			x$usedObj$pr = x$usedObj$prGenes = NULL
			invisible(x)
		} )
