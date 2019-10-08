context( 'pseudotime')
#library(BioData)
set.seed(1)
dat = matrix(round(rnorm(10000,mean = 1, sd = 1)),ncol=100)
dat[which(dat < 1)] = 0
colnames(dat) <- paste('Sample', 1:100)
rownames(dat) <- paste( 'gene', 1:100)
dat[,21:40] = dat[,21:40] + .2
dat[,41:60] = dat[,41:60] + .4
dat[,61:80] = dat[,61:80] + .6
dat[,81:100] = dat[,81:100] + .8

#dat[95:100,] = dat[95:100,] * -1

samples <- data.frame(SampleID = 1:100, sname = colnames(dat) )
annotation <- data.frame( GeneID = paste( 'gene', 1:100), Start= 101:200 )

x <- BioData$new( as.data.frame(dat), annotation=annotation,  Samples=samples, name="testObject",namecol='sname', outpath = "" )
x$samples$grouping = c( rep( 1, 20), rep(2, 20), rep( 3,20), rep( 4, 20), rep(5, 20))


x$outpath = tempdir()

mds( x, mds.type='UMAP', dim=2)
x$samples$grouping = c( rep( 1, 20), rep(2, 20), rep( 3,20), rep( 4, 20), rep(5, 20))
colors_4(x, 'grouping')
plotMDS( x, mds='Expression UMAP_2_dims', 'grouping', c() )

expect_true( file.exists( file.path(x$outpath, 'testObject_Expression UMAP_2_dims_1_2_grouping.pdf')) )

## And that does also look beautiful!
d= x$usedObj$MDS_PCA100_dim_2$'Expression UMAP_2_dims'
genes =pseudotimeTest( x, a = d[,1],b = d[,2], grouping = 'grouping', outpath = x$outpath, n = 5, smooth = 5, plotType = 'pdf')

expect_true(length(genes) == 10)

lapply( genes , function(n) {
			expect_true( file.exists( file.path(x$outpath, paste( n , sep="_","traj.corGenes_5_factor.pdf" )) ) )
})

unlink(file.path( x$outpath ,"*") )

mds( x, mds.type='UMAP', dim=3)

d= x$usedObj$MDS_PCA100$'Expression UMAP'

genes =pseudotimeTest3D( x, a = d[,1],b = d[,2],c= d[,3], grouping = 'grouping', outpath = x$outpath, n = 5, smooth = 5, plotType = 'pdf')




