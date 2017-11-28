#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @description create a NGSexpressionSet from a counts list object
#' @param dat the counts list you get from Rsubread::featureCounts()
#' @title description of function as_BioData
#' @export 
if ( ! isGeneric('as_BioData') ){ setGeneric('as_BioData', ## Name
	function ( dat ) { ## Argumente der generischen Funktion
		standardGeneric('as_BioData') ## der Aufruf von standardGeneric sorgt für das_BioData Dispatching
	}
)
}else {
	print ("Onload warn generic function 'as_BioData' already defined - no overloading here!")
}

setMethod('as_BioData', signature = c ('list'),
	definition = function ( dat ) {
	ret = NULL
	if (all.equal( names ( dat), c("counts" ,"annotation", "targets", "stat")  ) ) {
		samples <- data.frame(t(dat$stat))
		colnames(samples) <- as.character(t(samples[1,]))
		samples <- samples[-1,]
		samples$filename <- rownames(samples)
		rownames(samples) <-1:nrow(samples)
		ret <- BioData$new( 
				dat= cbind(dat$annotation, dat$counts), 
				Samples = samples, 
				namecol= 'filename', 
				namerow= 'GeneID',
				outpath= ''
		)
	}
	else {
		print ("The list needs to contain the entries counts ,annotation, targets and stat" )
	}
	ret
} )

setMethod('as_BioData', signature = c ('cellexalvr'),
		definition = function ( dat ) {
			dat <- cellexalvr::renew(dat)
			#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
			ok = which(lapply(colnames(dat$meta.cell) , function(x) { all.equal( as.character(as.vector(dat$meta.cell[,x])), colnames(dat$data)) == T } )==T)
			if ( length(ok) == 0) {
				if (all.equal( rownames(dat$meta.cell), colnames(dat$data)) ){
					dat$meta.cell = cbind( cell.name = colnames(dat$data), dat$meta.cell)
					namecol = 'cell.name'
				}
			}
			else {
				namecol = colnames(dat$meta.cell)[ok]
				namecol = namecol[1]
			}
			namerow = NULL
			if (nrow(dat$meta.gene)==0) {
				dat$meta.gene <- matrix(ncol=2, c(rownames(dat$data), rep( 0, nrow(dat$data)) ) )
				colnames(x@meta.gene) = c('Gene.Symbol', 'useless')
				rownames(x@meta.gene) = rownames(x@data)
				namerow = 'Gene.Symbol'
			}else {
				ok = which(lapply(colnames(dat$meta.gene) , function(x) { all.equal( as.character(as.vector(dat$meta.gene[,x])), rownames(dat$data)) == T } )==T)
				if ( length(ok) == 0) {
					if (all.equal( rownames(dat$meta.gene), rownames(dat$data)) ){
						dat$meta.gene = cbind( gene.name = rownames(dat$data), dat$meta.gene)
						namerow = 'gene.name'
					}
				}
				else {
					namerow = colnames(dat$meta.gene)[ok]
					namerow = make.names(namerow[1])
				}
			}
			storage.mode(dat$data) <- 'numeric'
			d <- data.frame(cbind( dat$meta.gene, dat$data))
			if ( nrow(dat$userGroups) == nrow(dat$meta.cell)){
				samples <- data.frame(cbind(dat$meta.cell, dat$userGroups))
			}else {
				samples <- data.frame(dat$meta.cell)
			}
			samples[,namecol] <- make.names(samples[,namecol])
			ret <- BioData$new( d, Samples=samples, name= 'from.cellexalvr', namecol= namecol, namerow=namerow, outpath='./' )
			ret$usedObj <- dat$usedObj
			ret
		}
		
)


setMethod('as_BioData', signature = c ('seurat'),
				definition = function ( dat ) {
					#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
					ok = which(lapply(colnames(dat@meta.data) , function(x) { all.equal( as.character(as.vector(dat@meta.data[,x])), colnames(dat@data)) == T } )==T)
					if ( length(ok) == 0) {
						if (all.equal( rownames(dat@meta.data), colnames(dat@data)) ){
							dat@meta.data = cbind( cell.name = colnames(dat@data), dat@meta.data)
							namecol = 'cell.name'
						}
					}
					else {
						namecol = colnames(dat@meta.data)[ok]
						namecol = namecol[1]
					}
					namerow = NULL
					if (nrow(dat@hvg.info)==0) {
						dat@hvg.info <- matrix(ncol=2, c(rownames(dat@data), rep( 0, nrow(dat@data)) ) )
						colnames(dat@hvg.info) = c('Gene.Symbol', 'useless')
						rownames(dat@hvg.info) = rownames(x@data)
						namerow = 'Gene.Symbol'
					}else {
						ok = which(lapply(colnames(dat@hvg.info) , function(x) { all.equal( as.character(as.vector(dat@hvg.info[,x])), rownames(dat@data)) == T } )==T)
						if ( length(ok) == 0) {
							if (all.equal( rownames(dat@hvg.info), rownames(dat@data)) == TRUE ){
								dat@hvg.info = cbind( gene.name = rownames(dat@data), dat@hvg.info)
								namerow = 'gene.name'
							}else {
								m <- match( rownames(dat@hvg.info), rownames(dat@data))
								dat@hvg.info = cbind( gene.name = rownames(dat@data), dat@hvg.info[order(m),])
								namerow = 'gene.name'
							}
						}
						else {
							namerow = colnames(dat@hvg.info)[ok]
							namerow = make.names(namerow[1])
							## now cast the hvg.info into the same order as the data object.
							m <- match( rownames(dat@hvg.info), rownames(dat@data))
							dat@hvg.info = cbind( gene.name = rownames(dat@data), dat@hvg.info[order(m),])
						}
					}
					#storage.mode(dat$data) <- 'numeric'
					d <- data.frame(cbind( dat@hvg.info, as.matrix(dat@data)))				
					samples <- data.frame(as.matrix(dat@meta.data))
					
					samples[,namecol] <- make.names(samples[,namecol])
					ret <- BioData$new( d, Samples=samples, name= 'from.cellexalvr', namecol= namecol, namerow=namerow, outpath='./' )
					ret$zscored <- as.data.frame(as.matrix(dat@scale.data))
					m <- match(colnames(dat@data), colnames(dat@raw.data))
					ret$raw <- as.data.frame(as.matrix(dat@raw.data))[,m]
					## now I need to manually remove unreliable data from the zscored information (set to -20)
					d <- function( i, obj) {
						drop <- which(obj$raw[,i] == 0)
						if ( length(drop) > 0 ){
							ret$zscored[drop,i] = -20
						}
						0
					}
					lapply(  1:ncol(ret$dat), d, ret)
					for ( i in 1:nrow(ret$dat) ) {
						drop <- which(ret$raw[i,] == 0)
						if ( length(drop) > 0 ){
							ret$zscored[i, drop] = -20
						}
					}
					
					#ret$usedObj <- dat$usedObj
					ret
				}  
)
