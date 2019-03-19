conext('Rcpp')

set.seed(1)
ncol=10
nrow=100
dat =  matrix(round(rnorm(nrow*ncol, 5, 10)),ncol=ncol)
bad = which(dat < 0 )
if ( length(bad) > 0 ) {
	dat[bad] = 0
}
bad = sample( 1:(ncol*nrow), 100 )
dat[bad] = -1

colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)


samples <- data.frame(SampleID = 1:ncol, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:nrow), Start= 101:(100+nrow) )

x <- BioData$new( cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = "" )

class(x) = c("SingleCells", "BioData", "R6")

z.score(x)

expect_equal( dim(x$dat), dim(x$zscored)) # dims are correct

zero  <- which( x$dat@x == -1 )
zeroZ <- which( x$zscored@x == -1 )

all_equal( zeroZ, zero) # zscore does not change the droped cells (-1)

d = x$zscored[1,]
expect_equal(mean(d[which(d > 0 )]), 10) # projects into 10+-1
expect_equal(sd(d[which(d > 0 )]), 1)

d = x$zscored[5,]
expect_equal(mean(d[which(d > 0 )]), 10) # projects into 10+-1
expect_equal(sd(d[which(d > 0 )]), 1)


