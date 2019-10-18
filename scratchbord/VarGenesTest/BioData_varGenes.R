library(BioData)
packageVersion('BioData')
#1] ‘1.3.5’

exp.mat <- read.table(file = "data/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
    as.is = TRUE, row.names = 1)

data = as_BioData( exp.mat )

class(data) = c('SingleCells', "BioData", "R6" )

data$samples$total=factor( rep(1, ncol(data$dat)) )
data$samples$nUMI = Matrix::colSums( data$dat)

genes = getGenesExpressedHigherThanExpected( data, 'total' )

cat( unlist(genes), file="BioData_varGenes_2000.txt" )


## once I am at it - why not check what my object can do?
normalize( data, 1e+5)
log.this(data)
z.score(data)

addCellCyclePhase( data )
CreateBin (data )

varGenes = SeuratVarGenes(data, 2000)

Intersect = names(which(table( c( unlist(genes), varGenes ) ) == 2) )
length(Intersect)
# 210

VarGenesObj <- reduceTo( data, what='row', to=Intersect, name='Intersect', copy=T)

mds( VarGenesObj, dim=2 )

plotMDS( VarGenesObj, mds="Raw Expression PCA_2_dims", c('Phase', 'nUMI_bins'), c('DNTT', 'ELANE', 'PROCR' ))

clusters( VarGenesObj , clusterby = "Raw Expression PCA_2_dims", groups.n = 10, ctype = "kmeans", onwhat= 'MDS', name = "FirstClustering" )

lapply( rev(c('FirstClustering', 'Phase', 'nUMI_bins')), function(n, obj){ reorder.samples(obj,n); colors_4(obj,n);NULL;}, VarGenesObj)
complexHeatmap( VarGenesObj, colGroups=c('FirstClustering', 'Phase', 'nUMI_bins'), ofile="FirstTest", pdf=T)

auto_reorder_grouping( VarGenesObj, 'FirstClustering')

complexHeatmap( VarGenesObj, colGroups=c('FirstClustering', 'Phase', 'nUMI_bins'), ofile="FirstTest_auto_reordered", pdf=T)

paste( collapse=", ",levels(VarGenesObj$samples$FirstClustering))
#[1] "5, 10, 6, 7, 2, 9, 1, 8, 3, 4"

new_order = list(
        "CD9_CD63_off" = c( 2, 9 ),
        "CD9_CD63_low" = c( 3, 4 ),
        "CD9_CD63_seconfHighest" = c( 1, 8 ),
        "CD9_CD63_high" = c(5, 10, 6, 7)
)

VarGenesObj$samples$FirstReorder = NULL
reorder_grouping( VarGenesObj, 'FirstClustering', unlist(new_order))
#[1] "sample group order changed"
#define_grouping( VarGenesObj,
complexHeatmap( VarGenesObj, colGroups=c('FirstClustering', 'Phase', 'nUMI_bins'), ofile="FirstTest_manual_reordered", pdf=T)



