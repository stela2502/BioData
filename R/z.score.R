#' @name z.score
#' @aliases z.score,BioData-method
#' @rdname z.score-methods
#' @docType methods
#' @description  z score the matrix
#' @param m the matrix of column = samples and rows = genes or a BioData object
#' @return the z scored matrix / BioData object
#' @title description of function z.score
#' @export 
if ( ! isGeneric('z.score') ){ methods::setGeneric('z.score', ## Name
		function (m) { 
			standardGeneric('z.score')
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
				m$zscored <- Matrix::Matrix(ret)
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
				m$zscored <- Matrix::Matrix(z.score( as.matrix( m$dat )))

			}
			invisible(m)
		})

setMethod('z.score', signature = c ('SingleCells'),
		definition = function ( m ) {
			
			if ( is.null(m$zscored) ) {
				m$zscored <- FastWilcoxTest::ZScore(m$dat)
				#ret[which(is.na(ret)==T)] <- 0
				colnames(m$zscored) <- colnames(m$dat)
				rownames(m$zscored) <- rownames(m$dat)
			}
			gc()
			invisible(m)
		})

