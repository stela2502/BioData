## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  cache=TRUE,
  autodep=TRUE, 
  comment = "#>",
  dataLoad = function(before, options, envir) {
  if (before == TRUE) {
    merged <- loadObj(file.path('TestData', 'BioData_resample_normalized_resample_normalized.RData') )
  } else if (before == FALSE)  {
    saveObj(merged)
  }
  },varLoad = function(before, options, envir) {
  if (before == TRUE) {
    merged <- loadObj(file.path('TestData', 'BioData_resample_normalized_resample_normalized.RData'))
    VarData <- loadObj(file.path('TestData', 'VariableDataset.RData'))
  } else if (before == FALSE)  {
    saveObj(merged)
    saveObj(VarData)
  }
  }, statLoad = function(before, options, envir) {
  if (before == TRUE) {
    merged <- loadObj(file.path('TestData', 'BioData_resample_normalized_resample_normalized.RData'))
    VarData <- loadObj(file.path('TestData', 'VariableDataset.RData'))
    StatGOI1 <- loadObj(file.path('TestData', 'SeuratStats1.RData'))
  } else if (before == FALSE)  {
    saveObj(merged)
    saveObj(VarData)
    saveObj(StatGOI1)
  }
  }
)

## ---- , echo=FALSE-------------------------------------------------------
library(BioData)

## ------------------------------------------------------------------------
m = matrix(0,ncol=10, nrow=10)
colnames(m) = paste('cell', 1:10)
rownames(m) = paste('gene', 1:10)
test <- as_BioData( m )
test
reduceTo(test, what='row', to=rownames(test$dat)[1:5])
test

## ------------------------------------------------------------------------
TM_Heart_10x <- as_BioData( system.file("extdata", "droplet_Heart_TabulaMuris.sqlite", package = "BioData"))
TM_Heart_10x$name = 'TM_Heart_10x'
TM_Heart_10x

## ---- merged, eval=FALSE-------------------------------------------------
#  TM_Aorta_facs <- as_BioData( system.file("extdata", "facs_Aorta_TabulaMuris.sqlite", package = "BioData"))
#  TM_Aorta_facs$name = 'TM_Aorta_facs'
#  TM_Aorta_facs
#  
#  merged <- merge( TM_Heart_10x, TM_Aorta_facs )
#  merged
#  

## ---- eval = FALSE-------------------------------------------------------
#  class(BioDataObject)

## ---- dataLoad="yes"-----------------------------------------------------
merged = TestData

merged$samples$expMethod =  
  unlist(lapply(str_split( merged$samples$samples, '\\.' ) , 
    function(x){ if (length(x) == 1) { 'dropseq' } else { 'smartseq' } }))

colors_4(merged, 'expMethod')

merged$outpath = file.path(getwd(), 'TestData')
if ( ! file.exists( merged$outpath) )
  dir.create( merged$outpath )
saveObj(merged)

hist(log10(TestData$samples$nUMI))

## ---- dataLoad="yes"-----------------------------------------------------
normalize(merged, reads=2500 )
z.score(merged)

normalized_reads = apply (merged$dat, 2, function(x) { OK = which( x > 0 ); if ( length(OK) > 0 ) { sum(x[OK])} else {0} }  )
all( normalized_reads == 2500)


## ---- dataLoad="yes"-----------------------------------------------------

addCellCyclePhase(merged, gnameCol='gname' )

table(merged$samples$Phase)

colors_4(merged, 'Phase')


## ---- dataLoad="yes"-----------------------------------------------------

V = apply( merged$dat, 1, function(x) { 
    OK = which(x != -1 ) 
    if ( length(OK) > 3) { var(x[OK])} 
    else {1e-5} } 
)

hist(log10(V))

V = V[rev(order(V))]
OK = names(V[1:700])
length(OK)

# Get rid of genes known/likely to cause problems
lapply( c( 'Mt-', 'Rp', 'GM', 'Rik' ) , function(n) { 
  bad = grep(n, OK )
  if ( length(bad ) > 0 ) {
    OK = OK[-bad]
  }
  NULL;
} )

length(OK)

CreateBin( merged )

colors_4(merged, 'expMethod')

VarData =reduceTo(merged, copy=T, what='row', to=OK, name="VariableDataset" )
saveObj(VarData)
VarData


## ---- varLoad='yes'------------------------------------------------------

mds( VarData, mds.type='TSNE_R', useRaw=T )
clusters( VarData, 
  clusterby= 'TSNE_R', groups.n= 10, 
  ctype='kmeans', onwhat='MDS',  
  name="First_TSNE_grouping" )


## ---- firstHeatmap, varLoad='yes', figure.height = 1000, figure.width = 1000----

lapply( 
  rev(c('First_TSNE_grouping', 'expMethod', 'Phase', 'nUMI_bins') ), 
  function(n) { reorder.samples( VarData, n ); NULL; } 
  )
dev.off()
complexHeatmap( VarData, colGroups=c('First_TSNE_grouping', 'expMethod', 'Phase', 'nUMI_bins'))
dev.off()

## ---- secondHeatmap, varLoad='yes', figure.height = 1000, figure.width = 1000----

m = match( colnames(merged$dat), colnames(VarData$dat))
all.equal(  colnames(merged$dat), colnames(VarData$dat)[m])

merged$samples$First_TSNE_grouping = VarData$samples$First_TSNE_grouping[m]

Seurat_FindAllMarkers( merged, 'First_TSNE_grouping', logfc.threshold = 5  )

names(merged$stats$Seurat_FindAllMarkers_First_TSNE_grouping)

ret_genes = 30
top_genes <- function( x ) {
  if ( length(x) == 0) {
    NA
  }
  else if ( length(x) < ret_genes ) {
    x
  }else {
    x[1:ret_genes]
  }
}

genes_list <- split( 
  as.vector(merged$stats$Seurat_FindAllMarkers_First_TSNE_grouping[,'gene']), 
  merged$stats$Seurat_FindAllMarkers_First_TSNE_grouping[,'cluster'] 
  )
deg.genes = unique(unlist( lapply( genes_list,top_genes ) ))
colors_4(merged, 'First_TSNE_grouping')

StatGOI1 <- reduceTo( merged, copy= T, what='row', to=deg.genes, name="SeuratStats1" )

lapply( 
  rev(c('First_TSNE_grouping', 'expMethod', 'Phase', 'nUMI_bins') ), 
  function(n) { reorder.samples( StatGOI1, n ); NULL; } 
  )
dev.off()
complexHeatmap( StatGOI1, 
  colGroups=c('First_TSNE_grouping', 'expMethod', 'Phase', 'nUMI_bins'))
dev.off()
saveObj(StatGOI1)


## ---- statLoad='yes'-----------------------------------------------------

dev.off()
complexHeatmap( 
  StatGOI1, 
  colGroups=c('First_TSNE_grouping', 'expMethod', 'Phase', 'nUMI_bins'), 
  green=T)


## ---- StatGOI1,  statLoad='yes', figure.height = 1000, figure.width = 1000----

mds( StatGOI1, mds.type='TSNE_R', useRaw=T )
clusters( 
 StatGOI1, clusterby= 'TSNE_R', 
 groups.n= 10, ctype='kmeans', 
 onwhat='MDS',  name="Second_TSNE_grouping" )

dev.off()
complexHeatmap( 
  StatGOI1, 
  colGroups=c('Second_TSNE_grouping', 'expMethod', 'Phase', 'nUMI_bins'), 
  green=T)


## ---- statLoad='yes', figure.height = 1000, figure.width = 1000----------
transpose ( StatGOI1 )

mds( StatGOI1, mds.type='TSNE_R', useRaw=T )
clusters( 
  StatGOI1, 
  clusterby= 'TSNE_R', 
  groups.n= 10, 
  ctype='kmeans', 
  onwhat='MDS',  
  name="Second_TSNE_Gene_grouping" )

transpose ( StatGOI1 )

dev.off()
complexHeatmap( 
   StatGOI1,
   colGroups=c('Second_TSNE_grouping', 'expMethod', 'Phase', 'nUMI_bins'),
   rowGroups= c('Second_TSNE_Gene_grouping'), 
   green=T)


## ---- statLoad='yes'-----------------------------------------------------
lapply(levels( StatGOI1$annotation$Second_TSNE_Gene_grouping ),
  function(n) {
    sort( rownames(StatGOI1$dat)[
      which( StatGOI1$annotation$Second_TSNE_Gene_grouping == n) ] 
     )  
  }
)

## ---- statLoad='yes'-----------------------------------------------------

lapply(levels( StatGOI1$annotation$Second_TSNE_Gene_grouping ), 
  function(n) { 
	  genes = sort( rownames(StatGOI1$dat)[ 
    which( StatGOI1$annotation$Second_TSNE_Gene_grouping == n) ] 
                    )
	  paste( sep="", 'http://david.abcc.ncifcrf.gov/api.jsp?type=','OFFICIAL_GENE_SYMBOL',"&ids=", paste(collapse=",", genes), "&tool=summary" )
        } 
)


## ------------------------------------------------------------------------
sessionInfo()

