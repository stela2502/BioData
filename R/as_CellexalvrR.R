#' @name as_cellexalvr
#' @aliases as_cellexalvr,BioData-method
#' @rdname as_cellexalvr-methods
#' @docType methods
#' @description export a BioData as cellexalVR folder that can be loaded into the CellexalVR VR platform
#' @param x The BioData object
#' @param outpath the outpath for the VR data
#' @param specie which specie the data comes from (mouse or human)
#' @param meta.cell.groups which BioData samples groups should become available in the Attributes coloring mode in VR
#' @param meta.genes.groups at the moment saver to leave NULL
#' @param userGroups  which sample columns should become accessible as userGroupings in VR to re-color the graphs
#' @param minGene VR has a problem with 0 level genes, hence these need to be removes, but you can also set a higher cut off (default=1)
#' @title description of function as_cellexalvr
#' @export 
setGeneric('as_cellexalvr', ## Name
	function (x, outpath, specie, meta.cell.groups, meta.genes.groups = NULL, userGroups=NULL, minGene=1 ) { ## Argumente der generischen Funktion
		standardGeneric('as_cellexalvr') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('as_cellexalvr', signature = c ('BioData'),
	definition = function (x, outpath, specie, meta.cell.groups, meta.genes.groups = NULL, userGroups=NULL, minGene=1 ) {
	normUMI <- apply( x$dat, 1, function(d) { sum(exp(d[which(d > 0 )] )) } )
	forCellexal <- reduceTo( x, what='row', to = rownames(merged$dat)[which(normUMI >= minGene  ) ], name=x$name )
	## get rid of -1 values
	bad = which(forCellexal$dat@x < 0 )
	if ( length(bad) > 0 )
		forCellexal$dat@x[bad] = 0
	
	## save space
	forCellexal$raw = NULL
	forCellexal$zscored = NULL
	forCellexal$dat = Matrix::drop0(forCellexal$dat)
	
	## store the important variables for the convert
	forCellexal$usedObj$meta.cell.groups = meta.cell.groups
	forCellexal$usedObj$specie = specie
	forCellexal$usedObj$meta.genes.groups = meta.genes.groups
	forCellexal$usedObj$userGroups = userGroups
	forCellexal$usedObj$outpath = outpath
	forCellexal$usedObj$specie = specie
	
	save( forCellexal, file="forCellexal.RData" )
	
	script = paste( sep="\n", "library(cellexalvrr)",
			"options(warn=-1)",
			paste( sep="", 'load("forCellexal.RData")' ),
			paste( sep="", "try({class(forCellexal) = 'list' }, silent=T) ## changes class to environment"  ),
			paste( sep="", 
"test = cellexalvrr( forCellexal, 
	meta.cell.groups =forCellexal$usedObj$meta.cell.groups,
	meta.genes.groups = forCellexal$usedObj$meta.genes.groups,
	userGroups = forCellexal$usedObj$userGroups,
	outpath = forCellexal$usedObj$outpath,
	specie = forCellexal$usedObj$specie )" ),
			paste( sep="", "if ( !file.exists(forCellexal$usedObj$outpath)) { dir.create(forCellexal$usedObj$outpath)}" ),
			paste( sep="", "export2cellexalvr( test, '", forCellexal$usedObj$outpath ,"' )")
	)
	print (paste ("create and run script forCellexal_convert.R" ) )
	system( "Rscript forCellexal_convert.R")
	print ( paste("please check the outpath", forCellexal$usedObj$outpath ))
	invisible(x)
} )
