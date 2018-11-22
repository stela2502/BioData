#library(BioData)
context( "single cells" )

set.seed(1)
norm <- matrix(rnorm(1000, mean=20),ncol=10)
null_indx <- c(
		62, 355, 576, 534, 602, 484, 131, 958, 737, 364, 
		150, 225, 992, 50, 586, 151, 751, 979, 171, 608, 
		528, 802, 315, 569, 70, 557, 547, 423, 936, 691, 
		251, 897, 408, 790, 730, 398, 496, 523, 938, 71, 
		739, 837, 949, 924, 609, 63, 123, 470, 632, 181, 
		169, 110, 149, 812, 402, 831, 779, 957, 598, 651, 
		745, 33, 770, 645, 405, 229, 8, 393, 734, 555, 917, 
		491, 486, 966, 82, 237, 330, 230, 162, 980, 445, 
		940, 689, 319, 740, 231, 542, 143, 160, 207, 721, 
		354, 828, 111, 895, 684, 281, 18, 616, 396
)
norm[null_indx] = 0
norm[null_indx[1:3]] = -1 # lost in the norm process
dat = Matrix( norm )


colnames(dat) <- paste('Sample', 1:10)
rownames(dat) <- paste( 'gene', 1:100)
samples <- data.frame(SampleID = 1:10, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )
dat <- round(dat)

x <- SingleCells$new( dat, annotation=annotation, Samples=samples, name="testObject",namecol='sname', outpath = "" )

expect_equal( class( x ), c( "SingleCells", 'BioData', 'R6') )

## test the renew
x <- renew(x)

expect_equal( class( x ), c( "SingleCells", 'BioData', 'R6') , info="class OK after renew")

expect_equal( x$data(),x$dat, info ="the data accessory function is OK" ) 
## z.score 

context( "single cells z.score" )

z.score(x)

expect_true ( length(which(as.matrix(x$zscored)[null_indx[4:length(null_indx)]] != 0)) == 0, "null normalized as 0")

expect_true ( length(which(as.matrix(x$zscored)[null_indx[1:3]] != -1)) == 0, "-1 values keept")

expect_equal( x$data(), x$zscored, info = "the data accessory function returns the z.scored data" ) 

expect_true ( round(sd(x$zscored@x[which( x$zscored@x > 0)])) == 1, paste("z.scored sd not 1!", sd(x$zscored@x[which( x$zscored@x > 0)])) )

context( "single cells mds" )

mds(x)

expect_true( ! is.null(x$usedObj$pr), "initial pca speed up not existing" )

expect_true( ! is.null(x$usedObj$MDS_PCA100[['Expression PCA']]), "PCA data not existing" )

