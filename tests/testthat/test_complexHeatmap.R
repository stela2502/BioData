test_that("heatmap",{
set.seed(1)
dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
colnames(dat) <- paste('Sample', 1:10)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:10, sname = colnames(dat), Cgroup= rep(c(1,2),5) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200, Rgroup= rep(1:5,20) )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "testObject" )
colors_4(x, c( 'Cgroup', 'Rgroup'))

complexHeatmap( x, colGroup=c('Cgroup') , ofile="Test1" )


complexHeatmap( x, colGroup=c('Cgroup') , rowGroup=c('Rgroup'),ofile="Test2", pdf=T )

complexHeatmap( x, colGroup=c('Cgroup') ,ofile="Test3", pdf=T )

complexHeatmap( x, ofile="Test4", pdf=T )

files <- c( 'testObject/Test1.png', 'testObject/testObject_col_Cgroup.pdf', 'testObject/testObject._legend_values.pdf',
'testObject/Test2.pdf','testObject/Test3.pdf', 'testObject/Test4.pdf', 'testObject/testObject_col_Cgroup.png', 'testObject/testObject_row_Rgroup.pdf' )
for ( f in  files ) {
	expect_true( file.exists(f), f)
	if (file.exists(f) ) { file.remove( f ) }
}

})