context( 'mds')
set.seed(1)
dat = matrix(round(rnorm(10100,mean = 10, sd = 15)),ncol=101)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:101)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:101, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )

## now I have a 10x101 data table which I can mds and check the outcome.
context( "mds PCA" )
			
	mds( x )
	expect_equal (names(x$usedObj$MDS_PCA100), c('Expression PCA'))
	expect_equal (dim(x$usedObj$MDS_PCA100$'Expression PCA'), c(101,  3)) # 10 samples
	
	mds( x , genes=T )
	expect_equal (names(x$usedObj$MDSgene), c('Expression PCA'))
	expect_equal (dim(x$usedObj$MDSgene$'Expression PCA'), c(100,  3)) # 10 samples
	
skip ("skip DM LLE and DDRTree - too simple data?!")

context( "mds DM" )

			mds( x , mds.type='DM')
			expect_equal (names(x$usedObj$MDS_PCA100), c('Expression PCA', 'Expression DM'))
			expect_equal (dim(x$usedObj$MDS_PCA100$'Expression DM'), c(101,  3)) # 10 samples
			
			mds( x , genes=T , mds.type='DM')
			expect_equal (names(x$usedObj$MDSgene), c('Expression PCA', 'Expression DM'))
			expect_equal (dim(x$usedObj$MDSgene$'Expression DM'), c(100,  3)) # 10 samples
				

context( "mds LLE" )
			
			mds( x , mds.type='LLE', LLEK=5)
			expect_equal (names(x$usedObj$MDS_PCA100), c('Expression PCA', 'Expression DM', 'Expression LLE'))
			expect_equal (dim(x$usedObj$MDS_PCA100$'Expression LLE'), c(101,  3)) # 10 samples
			
			mds( x , genes=T , mds.type='LLE', LLEK=5)
			expect_equal (names(x$usedObj$MDSgene), c('Expression PCA', 'Expression DM', 'Expression LLE'))
			expect_equal (dim(x$usedObj$MDSgene$'Expression LLE'), c(100,  3)) # 10 samples
			

# not working
#test_that( "mds ZIFA" ,{
#			
#			mds( x , mds.type='ZIFA')
#			expect_equal (names(x$usedObj$MDS), c('Expression PCA', 'Expression DM', 'Expression LLE', 'Expression ZIFA'))
#			expect_equal (dim(x$usedObj$MDS$'Expression ZIFA'), c(10,  3)) # 10 samples
#			
#			mds( x , genes=T , mds.type='ZIFA')
#			expect_equal (names(x$usedObj$MDSgene), c('Expression PCA', 'Expression DM', 'Expression LLE', 'Expression ZIFA'))
#			expect_equal (dim(x$usedObj$MDSgene$'Expression ZIFA'), c(101,  3)) # 10 samples
#			
#		} )## "mds ZIFA"


context( "mds DDRTree")
			
			mds( x , mds.type='DDRTree')
			expect_equal (names(x$usedObj$MDS_PCA100), c('Expression PCA', 'Expression DM', 'Expression LLE', 'Expression DDRTree'))
			expect_equal (dim(x$usedObj$MDS_PCA100$'Expression DDRTree'), c(101,  3)) # 10 samples
			
			mds( x , genes=T , mds.type='DDRTree')
			expect_equal (names(x$usedObj$MDSgene), c('Expression PCA', 'Expression DM', 'Expression LLE', 'Expression DDRTree'))
			expect_equal (dim(x$usedObj$MDSgene$'Expression DDRTree'), c(100,  3)) # 10 samples
			
