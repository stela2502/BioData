
#library(BioData)
context("collaps")
set.seed(1)
dat = data.frame( matrix(rnorm(1000),ncol=10) ) 
colnames(dat) <- paste('Sample', 1:10)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )
x$samples$group <- rep(c('A','B'), 5)

normalize(x, to=apply( x$dat,2, sum))
z.score(x)

a <- collaps( x, groupCol='group', by='median')

expect_equal( dim(a$data()), c(100,2),info= "right drop" )

expect_equal( colnames(a$data()), c('A','B'), info="colnames")

## I am not going to check the data here - too much work for now.
