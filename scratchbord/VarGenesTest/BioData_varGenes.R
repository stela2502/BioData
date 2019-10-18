library(BioData)
packageVersion('BioData')
#1] ‘1.3.5’

exp.mat <- read.table(file = "data/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
    as.is = TRUE, row.names = 1)

data = as_BioData( exp.mat )

class(data) = c('SingleCells', "BioData", "R6" )

data$samples$total=factor( rep(1, ncol(data$dat)) )
data$samples$nUMI = Matrix::colSums( data$dat)

genes = getGenesExpressedHigherThanExpected( data, 'total' )

cat( unlist(genes), file="BioData_varGenes_2000.txt" )

