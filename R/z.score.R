#' @name z.score
#' @aliases z.score,BioData-method
#' @rdname z.score-methods
#' @docType methods
#' @description  z score the matrix
#' @param m the matrix of column = samples and rows = genes or a BioData object
#' @return the z scored matrix / BioData object
#' @title description of function z.score
#' @export 
setGeneric('z.score', ## Name
		function (m) { ## Argumente der generischen Funktion
			standardGeneric('z.score') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)


setMethod('z.score', signature = c ('tRNAMINT'),
		definition = function ( m ) {
			
			if ( ! m$zscored ) {
				m$raw <- m$data
				ma  <- as.matrix(m$data)
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
				m$data <- data.frame(ret)
				colnames(m$data)<- colnames(m$raw)
				m$zscored = TRUE
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
			if (! m$zscored ){
				m$raw <- m$data
				m$data <- data.frame(z.score( as.matrix( m$data )))
				m$zscored = TRUE
			}
			invisible(m)
		})
