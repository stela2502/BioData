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
setGeneric('reduceTo', ## Name
		function ( x, what='row', to=NULL, ... ) { ## Argumente der generischen Funktion
			standardGeneric('reduceTo') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('reduceTo', signature = c ('BioData'),
		definition = function ( x, what='row', to=NULL, name='reducedTo' ) {
			if ( ! is.null(to)) {
				if ( what =="row") {
					if ( length(which(is.na(match(to,rownames(x$data)))==F)) > 0 ) {
						useOnly <- match(to, rownames(x$data))
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
						x$data <- x$data[useOnly,]
						x$annotation <- x$annotation[useOnly,]
						
						if ( ! is.null( x$raw) ) {
							x$raw <- x$raw[useOnly,]
						}
						if ( length(x$stats) > 0 ) {
							lapply( names(x$stats) , function(name) {
										x$stats[[name]] <- x$stats[[name]][match(tolower(to) ,tolower(x$stats[[name]][,1]) ),]
									} )
						}
						x$name = name
					}else {
						print (paste( "None of the probesets matched the probesets in",x$objects[[name]]@name, "-> keep everything!"))
					}
					
					
				}else if ( what =="col" ) {
					to <- tolower(make.names(to))
					if ( length(which(is.na(match(to,tolower(colnames(x$data))))==F)) > 0 ) {
						useOnly <- match(to, tolower(colnames(x$data)))
						not.matched <- to[is.na(useOnly)]
						if ( length(not.matched) > 0 ){
							print (paste('Problematic samples:', paste(not.matched,sep=', ')))
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
						x$data <- x$data[,useOnly]
						x$samples <- x$samples[useOnly,]
						
						if ( ! is.null( x$raw) ) {
							x$raw <- x$raw[,useOnly]
						}
						if ( length(x$stats) > 0 ) {
							x$stats = list()
						}
						x$name = name
						
					}else {
						print (paste( "None of the names (to) matched the sample names in",x$name, "-> keep everything!"))
					}
				}else {
					stop(paste( "the option what='",what,"' is not supported!", sep="") )
				}
			}
			invisible(x)
		} )