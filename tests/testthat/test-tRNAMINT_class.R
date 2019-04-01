context('tRNA class')

#library(BioData)
skip( "tRNA support depricated")
x <- loadObj('TestData.tRNAMINT.RData')

x <- renew(x)

x$outpath = "./"

expect_equal( dim(x$data), c(5000 ,24) )

extractCodonInformation(x)

exp <- c("X.MINTbase.Unique.ID", "tRF.sequence", "tRF.type.s.", 
		"Sequence.locations.in.tRNA.space..comma.deliminated.", 
		"reliability", "AlaAGC", "AlaCGC", "AlaTGC", "ArgACG", 
		"ArgCCG", "ArgCCT", "ArgTCG", "ArgTCT", "AsnATT", "AsnGTT", 
		"AspGTC", "CysACA", "CysGCA", "GlnCTG", "GlnTTG", "GluCTC", 
		"GluTTC", "GlyCCC", "GlyGCC", "GlyTCC", "HisGTG", "IleAAT", 
		"IleGAT", "IleTAT", "LeuAAG", "LeuCAA", "LeuCAG", "LeuTAA", 
		"LeuTAG", "LysCTT", "LysTTT", "MetCAT", "PheGAA", "ProAGG", 
		"ProCGG", "ProTGG", "SerACT", "SerAGA", "SerCGA", "SerGCT", 
		"SerTGA", "SupTTA", "ThrAGT", "ThrCGT", "ThrTGT", "TrpCCA", 
		"TrpTCA", "TyrGTA", "ValAAC", "ValCAC", "ValTAC")

expect_equal(colnames(x$annotation), exp )
expect_equal(x$usedObj$Codons, exp[6:length(exp)] )

res <- c( 8, 7, 8, 0, 0, 2, 0, 0, 1, 3, 9, 0, 4, 5, 3, 21, 15, 13, 19, 8, 3, 3, 2, 0, 2, 3, 5, 5, 2, 3, 3, 3, 2, 1, 1, 1, 2, 2, 1, 5, 4, 0, 1, 0, 0, 0, 0, 2, 5, 11, 14)
names(res) <- exp[6:length(exp)] 
value <- sampleCodonUsage(x, 'H9.PUS7KO.1_S6.all.reads' )

expect_equal(value, res)

plotCodonUsage(x,'H9.PUS7KO.1_S6.all.reads', fname="test" )

expect_equal(file.exists('test_pie.png'), TRUE)
if ( file.exists('test_pie.png') ) {
	unlink('test_pie.png')
}

expect_equal(file.exists('test_bars.png'), TRUE)
if ( file.exists('test_bars.png') ) {
	unlink('test_bars.png')
}





