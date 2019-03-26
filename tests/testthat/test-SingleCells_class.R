context( "single cells" )

set.seed(1)
ncol = 100
nrow = 90

norm <- matrix(round(rnorm(ncol * nrow, mean=20, sd= 30 )),ncol=ncol)
norm[which(norm < 0)] = 0

dat = Matrix::Matrix( norm, sparse=T )

colnames(dat) <- paste('Sample', 1:ncol)
rownames(dat) <- paste( 'gene', 1:nrow)
samples <- data.frame(SampleID = 1:ncol, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:nrow), Start= 101:101+nrow )


x <- SingleCells$new( dat, annotation=annotation, Samples=samples, name="testObject",namecol='sname', outpath = "" )

expect_equal( class( x ), c( "SingleCells", 'BioData', 'R6') )

## test the renew
x <- renew(x)



expect_equal( class( x ), c( "SingleCells", 'BioData', 'R6') , info="class OK after renew")

expect_equal( x$data(),x$dat, info ="the data accessory function is OK" ) 

context( "single cells z.score" )
## normalize 
reads = min(Matrix::colSums( x$dat))

old = x$clone()
normalize(x, reads)
#all_equal(as.character(dim(x$dat)) ,as.character(c(nrow,ncol)))
all_equal( as.vector(apply( x$dat,2, function(x) sum(x[which(x > 0 )]) )) , rep(reads, ncol(x$dat)) )

d= lapply ( 1:ncol(x$dat), function( i ) {
			raw = as.vector(x$raw[,i])
			norm = as.vector(x$dat[,i])
			lost = which(norm == -1 )
			norm_double = raw / sum(raw) * 100
			## have all values that have been lost been below 0 after the norm?
			if ( length(lost)>0 ){
				expect_equal( norm_double[lost] , rep( 0, length(lost)),1)
				## are all other values in the range +-1 of the norm data?
				expect_equal(norm[-lost], norm_double[-lost], 1 )
			}else {
				expect_equal(norm, norm_double, 1 )
			}
			NULL;
		} )


## z.score 
z.score(x)

expect_true ( length(which(as.matrix(x$zscored)[null_indx[4:length(null_indx)]] != 0)) == 0, "null normalized as 0")

expect_true ( length(which(as.matrix(x$zscored)[null_indx[1:3]] != -1)) == 0, "-1 values keept")

expect_equal( x$data(), x$zscored, info = "the data accessory function returns the z.scored data" ) 

expect_true ( round(sd(x$zscored@x[which( x$zscored@x > 0)])) == 1, paste("z.scored sd not 1!", sd(x$zscored@x[which( x$zscored@x > 0)])) )

context( "single cells mds" )

mds(x)

expect_true( is.null(x$usedObj$pr), "initial pca speed up existing" )

expect_true( is.null(x$usedObj$MDS[['Expression PCA']]), "PCA data existing" )

