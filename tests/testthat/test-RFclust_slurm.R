context( "SLURM RFclust" )
if ( ! system('sbatch --help 2>/dev/null ')==0 ) {
	skip( "no SLURM based system" )
}
set.seed(1)
dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
colnames(dat) <- paste('Sample', 1:10)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )

rfCluster_row( x, subset=79, rep=1, slice=1,k=3, nforest=5, ntree=5, settings = list( 'A' = "lsens2018-3-3", t="00:20:00" ,p="dell" ) )
rfCluster_col( x, subset=10, rep=1, slice=1,k=3, nforest=5, ntree=5, settings = list( 'A' = "lsens2018-3-3", t="00:20:00" ,p="dell" ) )

Sys.sleep( 15 )

rfCluster_row( x, subset=79, rep=1, slice=1,k=3, nforest=5, ntree=5, settings = list( 'A' = "lsens2018-3-3", t="00:20:00" ,p="dell" ) )
rfCluster_col( x, subset=10, rep=1, slice=1,k=3, nforest=5, ntree=5, settings = list( 'A' = "lsens2018-3-3", t="00:20:00" ,p="dell" ) )

expect_equal (levels(x$samples$'RFgrouping RFclust 2 n=3'), as.character(1:3) ) # 10 samples
expect_equal (levels(x$annotation$'RFgrouping RFclust_row 3'), as.character(1:3) ) # 10 samples