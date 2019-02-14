context( 'Class usage')
#library(BioData)
test_that("basic class usage",{
set.seed(1)
dat = matrix(round(rnorm(10000,mean = 1, sd = 1)),ncol=100)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:100)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:100, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )


x <- BioData$new( as.data.frame(dat), annotation=annotation,  Samples=samples, name="testObject",namecol='sname', outpath = "" )


context( 'mds TSNE_R')

mds( x, mds.type='TSNE_R')

expect_equal( class( x ), c('BioData', 'R6') )

context( 'reduceTo')

reduceTo(x,what='row', to=paste( 'gene', c(1,50:100)), name="rows dropped" )

expect_equal(  x$name, "rows dropped" )
expect_equal(dim(x$data()), c(52,100) )
expect_equal(dim(x$annotation), c(52,2))

expect_equal(rownames(x$data()), paste('gene', c(1,50:100)))

a <- x$clone()

reduceTo(x,what='col', to=paste( 'Sample', c(1,3,5,7), sep='.'), name="rows and cols dropped" )
expect_equal(  x$name, "rows and cols dropped" )
expect_equal(dim(x$data()), c(52,4) )
expect_equal(dim(x$annotation), c(52,2))
expect_equal(dim(x$samples), c(4,3))
expect_equal(colnames(x$data()), paste('Sample', c(1,3,5,7), sep='.'))

expect_equal(colnames(a$data()), paste('Sample', 1:100, sep='.'))
expect_equal(  a$name, "rows dropped" )

x$samples$group <- paste( 'Group', rep(c(1,2,3,4)) )


context( 'colors_4')

colors_4( x, 'group' )
exp = rainbow(4)
expect_equal( x$usedObj$colorRange$group, exp )
})


#x2 <- tRNAMINT$new( dat, annotation=annotation, Samples=samples, name="testObject2",namecol='sname', outpath = "" )
#expect_equal( class( x2 ), c('tRNAMINT', 'BioData', 'R6') )
#expect_equal(dim(x2$data()), c(100,10) ) ## got the original data



