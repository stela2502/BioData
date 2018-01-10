library(BioData)
set.seed(1)
dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
colnames(dat) <- paste('Sample', 1:10)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )

rfCluster_row( x, subset=79, rep=1, slice=1,k=3, nforest=50, ntree=50, settings = list( 'A' = "lsens2017-3-2", t="00:20:00" ,p="dell" ) )

Sys.sleep( 5 )

rfCluster_col( x, subset=10, rep=1, slice=1,k=3, nforest=50, ntree=50, settings = list( 'A' = "lsens2017-3-2", t="00:20:00" ,p="dell" ) )
