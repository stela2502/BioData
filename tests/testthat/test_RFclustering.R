library(BioData)
set.seed(1)
dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
colnames(dat) <- paste('Sample', 1:10)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )

rfCluster_row( x, subset=79, rep=1, slice=1,k=3, nforest=50, ntree=50 )
## this should not return but run the analysis - nevertheless I need to restart it once
rfCluster_row( x, subset=79, rep=1, slice=1,k=3, nforest=50, ntree=50 )

expect_equal(colnames(x$annotation), c('GeneID', 'Start', 'RFgrouping RFclust_row 1') )
expect_equal( names(x$usedObj$rfExpressionSets_row), c('testObject_RFclust_row_1'))
expect_equal( names(x$usedObj$rfObj_row), c('testObject_RFclust_row_1'))

z.score(x)

rfCluster_row( x, subset=79, rep=1, slice=1,k=3, nforest=50, ntree=50 , name="try2" )
Sys.sleep( 5 )
rfCluster_row( x, subset=79, rep=1, slice=1,k=3, nforest=50, ntree=50 , name="try2" )



expect_equal(colnames(x$annotation), c('GeneID', 'Start', 'RFgrouping RFclust_row 1', "RFgrouping try2 1") )


rfCluster_col( x, subset=10, rep=1, slice=1,k=3, nforest=50, ntree=50 )
Sys.sleep( 5 )
rfCluster_col( x, subset=10, rep=1, slice=1,k=3, nforest=50, ntree=50 )
	

expect_equal(colnames(x$samples), c('SampleID', 'sname', 'samples', 'RFgrouping RFclust 1' ) )