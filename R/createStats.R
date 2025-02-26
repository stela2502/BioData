#' @name createStats
#' @aliases createStats,NGSexpressionSet-method
#' @rdname createStats-methods
#' @docType methods
#' @description  calculate staistics on all possible groupings using the DEseq nbinomTest test Both
#' @description  together create the group "A vs. B"
#' @param x the NGSexpressionSet
#' @param condition the grouping column in the samples data.frame
#' @param files write the statistics tables (FALSE)
#' @param A the first component to analyze (Group A)
#' @param B the second component to analyze (Group B)
#' @return the NGSexpressionSet with a set of ststs tables
#' @title description of function createStats
#' @export 
if ( ! isGeneric('createStats') ){ methods::setGeneric('createStats', ## Name
	function (x, condition, files=F, A=NULL, B=NULL, ...) { 
		standardGeneric('createStats')
	}
)
}else {
	print ("Onload warn generic function 'createStats' already defined - no overloading here!")
}
#
setMethod('createStats', signature = c ('BioData'),
	definition = function (x, condition, files=F, A=NULL, B=NULL) {
		#stop( "Not implemented / broken!")
			if (!requireNamespace("DESeq2", quietly = TRUE)) {
				stop("DESeq2 needed for this function to work. Please install it.",
						call. = FALSE)
			}
		if ( nrow(x$data()) < 1e+3 ) {
			stop ( "Please calculate the statistics only for the whole dataset!" )
		}
		if ( is.na( match ( condition, colnames(x$samples))) ) {
			stop ( 'Please select a condition from the sample colnames' )
		}
		id <- match ( condition, names(x$usedObj[['cds']]))
		if ( is.na(id)) { id <- 1 } 
		x<- preprocess ( x, condition )
		conditions <- as.vector(unique(x$samples[,condition]))
		x$usedObj[['cds']][[id]] <- DESeq2::DESeq(x$usedObj[['cds']][[id]]) ## run it
		if ( ! is.null(A) && ! is.null(B)) {
			x <- add_to_stat ( x, 
				stat = as.data.frame(DESeq2::results(x$usedObj[['cds']][[id]], contrast= c( condition,  A, B ) )), 
				name = paste( A, B ,sep=' vs. ')
			)
		}
		else {
			for ( i in 1:(length(conditions)-1) ){
				for ( a in (i+1):length(conditions) ){
					x <- add_to_stat ( x, 
						stat = as.data.frame(
							DESeq2::results(x$usedObj[['cds']][[id]], 
							contrast= c( condition, conditions[i] , conditions[a] ) )
						), 
						name = paste( conditions[i], conditions[a],sep=' vs. ')
					)

				}
			}
		}
		if ( files ) {
			#writeStatTables( x )
		}
		#eval( detach( 'package:DESeq2' ) )
		x
})



setMethod('createStats', signature = c ( 'MicroArray') ,
		definition = function ( x, condition, files=F, A=NULL, B=NULL, form=NULL  ) {
			if (!requireNamespace("limma", quietly = TRUE)) {
				stop("limma needed for this function to work. Please install it.",
						call. = FALSE)
			}
			if ( ! class(x$samples[, condition]) == 'factor') {
				x$samples[, condition] <- as.factor( x$samples[, condition] )
			}
			str = levels(x$samples[,condition])
			Factor=x$samples[,condition]
			design <- stats::model.matrix( ~0+ Factor )
			colnames(design) <- stringr::str_replace_all( colnames(design), 'Factor', '' )
			fit <- limma::lm.series(data$data(), design)
			contr <- list()
			i=1
			str = levels(x$samples[,condition])
			cmps <- NULL
			for ( a in str[1:(length(str)-1)] ){
				for ( b in str[2:length(str)] ){
					contr[[paste( a,b,sep="-")]] = paste( a,b,sep="-") 
					cmps <- c( cmps, paste( a,b,sep="-") )
					i = i + 1
				}
			}
			contr$levels =design
			if ( !is.null(A) & ! is.null(B)) {
				contr = list()
				cmps <- NULL
				for ( a in A ) {
					for (b in B ) {
						contr[[paste(sep="-",A,B)]] = paste(sep="-",A,B)
						cmps = c(cmps, paste(sep="-",A,B))
					}
				}
				contr$levels =design
			}
			cont.matrix <- do.call(limma::makeContrasts, contr )
			fit2 <- limma::contrasts.fit(fit,cont.matrix)
			if ( is.null(x$stats) ) {
				x$stats <- list()
			}
			
			for ( i in cmps ) {
				x$stats[[i]] <- limma::topTable(fit2, coef = i, adjust='fdr',number=1000000)
				x$stats[[i]] <- x$stats[[i]][match(rownames(x$data()),rownames(x$stats[[i]])) ,]
				x$stats[[i]] <- cbind( rownames(x$data()), x$stats[[i]] )
			}
			#detach( 'package:limma' )
			invisible(x)
		})


setMethod('createStats', signature = c ( 'SingleCells') ,
		definition = function ( x, condition, files=F, A=NULL, B=NULL, covariates=NULL, form=NULL ) {
			
			stop("This package is depricated - please shift to using an other single cell package.",call. = FALSE)
			
		}
)


