#' @name z.score
#' @aliases z.score,BioData-method
#' @rdname z.score-methods
#' @docType methods
#' @description  z score the matrix
#' @param m the matrix of column = samples and rows = genes or a BioData object
#' @return the z scored matrix / BioData object
#' @title description of function z.score
#' @export 
if ( ! isGeneric('z.score') ){ setGeneric('z.score', ## Name
		function (m) { ## Argumente der generischen Funktion
			standardGeneric('z.score') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
}else {
	print ("Onload warn generic function 'z.score' already defined - no overloading here!")
}


setMethod('z.score', signature = c ('tRNAMINT'),
		definition = function ( m ) {
			
			if ( is.null(m$zscored) ) {
				#m$raw <- m$data
				ma  <- as.matrix(m$dat)
				i = 0
				opts <- unique(as.character(m$samples$NormalizationMode))
				
				norm.name <- function( x ) {
					for( name in opts) {
						x[m$samples$NormalizationMode == name] = norm.z ( as.numeric(x[m$samples$NormalizationMode == name]))
					}
					x
				}
				norm.z <-  function (x) {
					i = i+1
					n <- which(is.na(x) == T)
					if ( length(x) - length(n) > 1 ){
						if (length(n) == 0 ){
							x <-  scale(as.vector(t(x)))
						}
						else {
							x[-n] <- scale(as.vector(t(x[-n])))
							x[n] <- NA
						}
						
					}
					else {
						x[] = NA
					}
					x
				}
				
				ret <- t(
						apply(ma,1, norm.name)
				)
				#ret[which(is.na(ret)==T)] <- -20
				m$zscored <- Matrix(ret)
				colnames(m$zscored)<- colnames(m$dat)
			}
			invisible(m)
		})

setMethod('z.score', signature = c ('matrix'),
		definition = function (m ) {
			rn <- rownames( m )
			me <- apply( m, 1, mean )
			sd <- apply( m, 1, sd )
			sd[which(sd==0)] <- 1e-8
			m <- (m - me) /sd
			rownames(m) <- rn
			invisible(m)
		})

setMethod('z.score',signature = c ('BioData'),
		definition = function (m) {
			if ( is.null( m$zscored ) ){
				#m$raw <- m$data
				m$zscored <- Matrix(z.score( as.matrix( m$dat )))

			}
			invisible(m)
		})

setMethod('z.score', signature = c ('SingleCells'),
		definition = function ( m ) {
			
			if ( is.null(m$zscored) ) {
				i = 0
				m$zscored <- Matrix(t(
						apply(m$dat,1, function (x) {
									i = i+1
									n <- which(x <= 0)
									dropped = which(x == -1)
									if ( length(x) - length(n) > 1 ){
										if (length(n) == 0 ){
											x <-  10 + scale(as.vector(t(x)))
										}
										else {
											x[-n] <- 10 + scale(as.vector(t(x[-n])))
								#			x[n] <- 0
								#			if ( length(dropped) > 0) {
								#				x[dropped] <- -1
								#			}
										}
										
									}
									else {
								#		x[] = 0
									}
									x}
						)
				))
				#ret[which(is.na(ret)==T)] <- 0
				if ( length( which(is.na(m$zscored)) ) >0 ) {
					m$zscored[which(is.na(m$zscored))] <- 0
				}
				colnames(m$zscored) <- colnames(m$dat)
				rownames(m$zscored) <- rownames(m$dat)
			}
			gc()
			invisible(m)
		})

