if (!library("devtools", quietly = TRUE,logical.return=TRUE )) {
	install.packages(c('devtools'),  repos='https://ftp.acc.umu.se/mirror/CRAN/')
	library(devtools)
}
source("https://bioconductor.org/biocLite.R")
biocLite()
install()
