

#library(BioData)
#test_that("z.score",{
set.seed(1)
dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
colnames(dat) <- paste('Sample', 1:10)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )
expect_equal( class( x ), c('BioData', 'R6'), "correct class name" )

normalize(x)

context( 'z.score')
z.score(x)

expect_true( length( which( round(apply(x$zscored,1,mean),digits=10) != 0 )) == 0 ,"all means == 0" )

expect_true( length( which( round(apply(x$zscored,1,sd),digits=10) != 1 )) == 0 ,"all sd's == 1" )

#})