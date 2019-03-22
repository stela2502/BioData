#' @name normalize
#' @aliases normalize,BioData-method
#' @docType methods
#' @description  normalize the expression data (sample wise) using DEseq2
#' @param x The BioData object (NGS expression data, not single cells)
#' @param readCounts The number of reads from each bam file or another value you want to normalize the data to
#' @param to_gene_length FALSE whether or not to normalize the data to gene length
#' @param geneLengthCol the column in the annotation data.frame to (in addition) normalize the genes to (e.g. trancript length)
#' @param name the new name of the object (deafule old name + normalized)
#' @param force replace old norm data (FALSE)
#' @return the normalized data set (original data stored in NGS$raw
#' @title normalize a BioData::R6 object
#' @export normalize

if ( ! isGeneric('normalize') ){ setGeneric('normalize', ## Name
	function ( object, ... , name=NULL) { 
		standardGeneric('normalize')
	}
)
}else {
	print ("Onload warn generic function 'normalize' already defined - no overloading here!")
}

setMethod('normalize', signature = c ('BioData'),
		definition = function (  object, readCounts=NULL, to_gene_length=FALSE, geneLengthCol='transcriptLength', force=FALSE ,name=NULL) {
			if ( ! object$snorm ){
				
				if ( is.null(object$raw) ){
					object$raw <- object$dat
				}else {
					object$dat <- object$raw ## regen if e.g. force is used ;-)
				}
				if ( is.null( readCounts ) ) {
					readCounts <- as.vector( DESeq2::estimateSizeFactorsForMatrix ( as.matrix(object$raw)) )
				}
				object$samples$SizeFactor <- readCounts
				object$dat =  data.frame(t(apply(object$dat,1, function(a) { a / readCounts } ) ))
				colnames(object$dat) = colnames(object$raw)
				rownames(object$dat) = rownames(object$raw)
				if (length(to_gene_length ) > 0 ){
					for ( i in 1:nrow(object$dat)) {
						object$dat[i,] <- object$dat[i,]/ (object$annotation[i,geneLengthCol ] / 1000)
					}
				}
				object$snorm=TRUE
				if(is.null(name)){
					name = paste( object$name ,'DESeq2_normalized' )
				}
				object$name = name
				
				if ( object$logged) {
					object$logged = FALSE
					logThis(object) ## regen log
				}
			}
			invisible(object)
		})


#' @describeIn normalize a SingleCells::BioData::R6 object using subsampling
#' @docType methods
#' @description normalize the expression data by subsampling as described in PMID 24531970
#' @param x The SingleCells::BioData::R6 object
#' @param reads the required read depth
#' @param name the name of the new object
#' @param  force re-normalize this object (default FALSE)
#' @return the normalized data set (original data stored in slot 'raw'
#' @title normalize a SingleCells::BioData::R6 object
#' @export normalize
setMethod('normalize', signature = c ('SingleCells'),
		definition = function (  object, reads=600, force=FALSE , name=NULL) {
			if ( is.null( object$usedObj$snorm) ) {
				object$usedObj$snorm = 0
			}
			reads <- round(reads)
			
			normF <- function(x, n) {
				x <- as.vector(x)
				ok1 <- which( x > 0 )
				d <- sample(rep ( 1:n, x) , reads, replace=T)
				
				t <- table(d)
				
				dropped <- setdiff( ok1, as.numeric(names(t)) )
				x[as.numeric(names(t))] <- as.numeric(t)
				
				x[ dropped ] <- -1
				x
			}
			
			if ( force | object$usedObj$snorm == 0 ) {
				if ( length( object$samples$nUMI ) == 0 ) {
					object$samples$nUMI <- apply( object$dat, 2, sum)
				}
				if(is.null(name)){
					name = paste( object$name ,'resample_normalized' )
				}
				object <- reduceTo( object, what="col", to=as.character(object$samples[which(object$samples$nUMI >= reads), object$sampleNamesCol ]) 
						, name=name )
				
				if ( is.null(object$raw) ){
					object$raw <- object$dat
				}
				## resample the data
				n <- nrow(object$raw)
				object$dat[] <- 0
				#pb <- progress_estimated(100)
				#steps = ceiling(ncol(object$raw)/100)
				#object$dat <- FastSparseRowScale( object$raw, TRUE, FALSE, reads, TRUE)
				object$dat <- Matrix(apply( object$raw,2, normF, n ))
				rownames(object$dat) <- rownames(object$raw)
			}
			else {
				print ("Data was already normalized - skipped")
			}
			object$usedObj$snorm = 1
			invisible(object)
		}
)




#' @describeIn normalize normalize a MicoArray::BioData::R6 object using quantile normalization
#' @docType methods
#' @description  constructor that has to be implemented for a generic BioData
#' This generic version was meant for array data and I have not had the need nor time to implement this part.
#' @param x the MicoArray::BioData::R6 object
#' @param to a numeric vector to normalize the samples to. Has to have the same length as are columns in the data table 
#' @title normalize a MicoArray::BioData::R6 object
#' @export normalize
setMethod('normalize', signature = c ('MicroArray') ,
	definition = function ( object , to=NULL, name=NULL) {
		object$zscored=NULL
		df_rank <- apply(object$dat,2,rank,ties.method="min")
		df_sorted <- data.frame(apply(object$data(), 2, sort))
		df_mean <- apply(df_sorted, 1, mean)
		
		index_to_mean <- function(my_index, my_mean){
			return(my_mean[my_index])
		}
		object$samples$norm_to = 'quantile'
		
		df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
		rownames(df_final) <- rownames(df)
		
		object$dat = Matrix( as.matrix(df_final) )
		
		invisible(object)
})

