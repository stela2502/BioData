#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @description create a NGSexpressionSet from a counts list object
#' @param dat the counts list you get from Rsubread::featureCounts()
#' @title description of function as_BioData
#' @export 
setGeneric('as_BioData', ## Name
	function ( dat ) { ## Argumente der generischen Funktion
		standardGeneric('as_BioData') ## der Aufruf von standardGeneric sorgt für das_BioData Dispatching
	}
)

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
			#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
			ok = which(lapply(colnames(dat@meta.cell) , function(x) { all.equal( as.character(as.vector(dat@meta.cell[,x])), colnames(dat@data)) == T } )==T)
			if ( length(ok) == 0) {
				if (all.equal( rownames(dat@meta.cell), colnames(dat@data)) ){
					dat@meta.cell = cbind( cell.name = colnames(dat@data), dat@meta.cell)
					namecol = 'cell.name'
				}
			}
			else {
				namecol = colnames(dat@meta.cell)[ok]
				namecol = namecol[1]
			}
			namerow = NULL
			if (nrow(dat@meta.gene)==0) {
				dat@meta.gene <- matrix(ncol=2, c(rownames(dat@data), rep( 0, nrow(dat@data)) ) )
				colnames(x@meta.gene) = c('Gene.Symbol', 'useless')
				rownames(x@meta.gene) = rownames(x@data)
				namerow = 'Gene.Symbol'
			}else {
				ok = which(lapply(colnames(dat@meta.gene) , function(x) { all.equal( as.character(as.vector(dat@meta.gene[,x])), rownames(dat@data)) == T } )==T)
				if ( length(ok) == 0) {
					if (all.equal( rownames(dat@meta.gene), rownames(dat@data)) ){
						dat@meta.gene = cbind( gene.name = rownames(dat@data), dat@meta.gene)
						namerow = 'gene.name'
					}
				}
				else {
					namerow = colnames(dat@meta.gene)[ok]
					namerow = make.names(namerow[1])
				}
			}
			storage.mode(dat@data) <- 'numeric'
			d <- data.frame(cbind( dat@meta.gene, dat@data))
			ret <- BioData$new( d, Samples=data.frame(cbind(dat@meta.cell, dat@userGroups)), name= 'from.cellexalvr', namecol= namecol, namerow=namerow, outpath='./' )
			ret$usedObj <- dat@usedObj
			ret
		}  )