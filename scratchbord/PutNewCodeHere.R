as_CellexalvrR <- function(x, outpath, specie, meta.cell.groups, meta.genes.groups = NULL, userGroups=NULL, name="BioData_as_cellexalvrR", minGene=1 ) {
	normUMI <- apply( x$dat, 1, function(d) { sum(exp(d[which(d > 0 )] )) } )
	forCellexal <- reduceTo( x, what='row', to = rownames(merged$dat)[which(normUMI >= minGene  ) ], name=name )
	## get rid of -1 values
	bad = which(forCellexal$dat@x < 0 )
	if ( length(bad) > 0 )
		forCellexal$dat@x[bad] = 0
	
	## save space
	forCellexal$raw = NULL
	forCellexal$zscored = NULL
	forCellexal$dat = drop0(forCellexal$dat)
	
	## store the important variables for the convert
	forCellexal$usedObj$meta.cell.groups = meta.cell.groups
	forCellexal$usedObj$specie = specie
	forCellexal$usedObj$meta.genes.groups = meta.genes.groups
	forCellexal$usedObj$userGroups = userGroups
	forCellexal$usedObj$outpath = outpath
	forCellexal$usedObj$specie = specie
	
	save( forCellexal, file="forCellexal.RData" )
	
	script = paste( sep="\n", "library(cellexalvrR)",
			"options(warn=-1)",
			paste( sep="", 'load("forCellexal.RData")' ),
			paste( sep="", "try({class(forCellexal) = 'list' }, silent=T) ## changes class to environment"  ),
			paste( sep="", 
"test = cellexalvrR( forCellexal, 
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
}