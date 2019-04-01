
context('merge')
## For simplicity and checkability I will use two extraemely simple obejcts
set.seed(1)
ncol= 1000
nrow=900

m1 = matrix ( 1, nrow=nrow, ncol=ncol)
colnames (m1) = paste("cell_1", sep=".", 1:ncol)
rownames(m1) = paste( "gene" , sep=".", 1:nrow)

d1 = as_BioData( m1 )


m2 = matrix ( 2, nrow=nrow, ncol=ncol)
colnames (m2) = paste("cell_2", sep=".", 1:ncol)

rownames(m2) = paste( "gene" , sep=".", sample( 1:(2*nrow), nrow ) )

d2 = as_BioData( m2 )

merged = merge( d1,d2)

expect_equal( as.vector(merged$dat[1,]), c(rep(1,ncol ), rep(0, ncol) ) )
expect_equal( as.vector(merged$dat[3,]), c(rep(1,ncol ), rep(2, ncol) ) )

expect_equal( dim(merged$dat) , c( 1341,2000) )
