dataframe2biodata <- function(x) {
	Samples <- data.frame( SampleName = colnames(x)[-1] )
	f <- str_replace(f, "\\.\\w\\w\\w",'')
	ret <- BioData$new( dat=x, Sample=Samples, namecol='SampleName', namerow=  colnames(x)[1], outpath=pwd(), name=f)
	ret
}


