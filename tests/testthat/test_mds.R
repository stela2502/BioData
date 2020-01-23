#library(BioData)

set.seed(1)
dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
colnames(dat) <- paste('Sample', 1:10)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )

## now I have a 10x100 data table which I can mds and check the outcome.
test_that( "mds PCA" ,{
			
	mds( x )
	expect_equal (names(x$usedObj$MDS), c('Expression PCA'))
	expect_equal (dim(x$usedObj$MDS$'Expression PCA'), c(10,  3)) # 10 samples
	
	mds( x , genes=T )
	expect_equal (names(x$usedObj$MDSgene), c('Expression PCA'))
	expect_equal (dim(x$usedObj$MDSgene$'Expression PCA'), c(100,  3)) # 10 samples
	
} )## "mds PCA"

test_that( "mds DM" ,{
			
			mds( x , mds.type='DM')
			expect_equal (names(x$usedObj$MDS), c('Expression PCA', 'Expression DM'))
			expect_equal (dim(x$usedObj$MDS$'Expression DM'), c(10,  3)) # 10 samples
			
			mds( x , genes=T , mds.type='DM')
			expect_equal (names(x$usedObj$MDSgene), c('Expression PCA', 'Expression DM'))
			expect_equal (dim(x$usedObj$MDSgene$'Expression DM'), c(100,  3)) # 10 samples
			
} )## "mds DM"

test_that( "mds LLE" ,{
			
			mds( x , mds.type='LLE', LLEK=5)
			expect_equal (names(x$usedObj$MDS), c('Expression PCA', 'Expression DM', 'Expression LLE'))
			expect_equal (dim(x$usedObj$MDS$'Expression LLE'), c(10,  3)) # 10 samples
			
			mds( x , genes=T , mds.type='LLE', LLEK=5)
			expect_equal (names(x$usedObj$MDSgene), c('Expression PCA', 'Expression DM', 'Expression LLE'))
			expect_equal (dim(x$usedObj$MDSgene$'Expression LLE'), c(100,  3)) # 10 samples
			
		} )## "mds LLE"

# not working
#test_that( "mds ZIFA" ,{
#			
#			mds( x , mds.type='ZIFA')
#			expect_equal (names(x$usedObj$MDS), c('Expression PCA', 'Expression DM', 'Expression LLE', 'Expression ZIFA'))
#			expect_equal (dim(x$usedObj$MDS$'Expression ZIFA'), c(10,  3)) # 10 samples
#			
#			mds( x , genes=T , mds.type='ZIFA')
#			expect_equal (names(x$usedObj$MDSgene), c('Expression PCA', 'Expression DM', 'Expression LLE', 'Expression ZIFA'))
#			expect_equal (dim(x$usedObj$MDSgene$'Expression ZIFA'), c(100,  3)) # 10 samples
#			
#		} )## "mds ZIFA"


test_that( "mds DDRTree" ,{
			
			mds( x , mds.type='DDRTree')
			expect_equal (names(x$usedObj$MDS), c('Expression PCA', 'Expression DM', 'Expression LLE', 'Expression DDRTree'))
			expect_equal (dim(x$usedObj$MDS$'Expression LLE'), c(10,  3)) # 10 samples
			
			mds( x , genes=T , mds.type='DDRTree')
			expect_equal (names(x$usedObj$MDSgene), c('Expression PCA', 'Expression DM', 'Expression LLE', 'Expression DDRTree'))
			expect_equal (dim(x$usedObj$MDSgene$'Expression DDRTree'), c(100,  3)) # 10 samples
			
		} )## "mds DDRTree"