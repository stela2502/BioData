#' @name createStats
#' @aliases createStats,NGSexpressionSet-method
#' @rdname createStats-methods
#' @docType methods
#' @description  calculate staistics on all possible groupings using the DEseq nbinomTest test Both
#' @description  together create the group 'A vs. B'
#' @param x the NGSexpressionSet
#' @param condition the grouping column in the samples data.frame
#' @param files write the statistics tables (FALSE)
#' @param A the first component to analyze (Group A)
#' @param B the second component to analyze (Group B)
#' @return the NGSexpressionSet with a set of ststs tables
#' @title description of function createStats
#' @export 
if ( ! isGeneric('createStats') ){ setGeneric('createStats', ## Name
	function (x, condition, files=F, A=NULL, B=NULL) { ## Argumente der generischen Funktion
		standardGeneric('createStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)
}else {
	print ("Onload warn generic function 'createStats' already defined - no overloading here!")
}
#
setMethod('createStats', signature = c ('BioData'),
	definition = function (x, condition, files=F, A=NULL, B=NULL) {
		#stop( "Not implemented / broken!")
		
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
						stat = as.data.frame(DESeq2::results(x$usedObj[['cds']][[id]], contrast= c( condition, conditions[i] , conditions[a] ) )), 
						name = paste( conditions[i], conditions[a],sep=' vs. ')
					)

				}
			}
		}
		if ( files ) {
			writeStatTables( x )
		}
		x
})

add_to_stat <- function( x, stat, name ) {
	if ( ! is.na( match( name, names(x$stats)))){
		x$stats[[ match( name, names(x$stats)) ]] <- stat
	}else {
		x$stats[[ length( x$stats ) +1 ]] <- stat
		names(x$stats)[length(x$stats) ] <- name
	}
	x
}


setMethod('createStats', signature = c ( 'MicroArray') ,
		definition = function ( x, condition, files=F, A=NULL, B=NULL ) {
			if (!requireNamespace("limma", quietly = TRUE)) {
				stop("limma needed for this function to work. Please install it.",
						call. = FALSE)
			}
			if ( ! class(x$samples[, condition]) == 'factor') {
				x$samples[, condition] <- as.factor( x$samples[, condition] )
			}
			str = levels(x$samples[,condition])
			Factor=x$samples[,condition]
			design <- model.matrix( ~ -1+ Factor )
			colnames(design) <- str_replace_all( colnames(design), 'Factor', '' )
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
			cont.matrix <- do.call(limma::makeContrasts, contr )
			fit2 <- limma::contrasts.fit(fit,cont.matrix)
			if ( is.null(x$stats) ) {
				x$stats <- list()
			}
			
			for ( i in cmps ) {
				x$stats[[i]] <- limma::toptable(fit2, coef = i, adjust='fdr',number=1000000)
				x$stats[[i]] <- x$stats[[i]][match(rownames(x$data()),rownames(x$stats[[i]])) ,]
				x$stats[[i]] <- cbind( rownames(x$data()), x$stats[[i]] )
			}
			invisible(x)
		})


setMethod('createStats', signature = c ( 'SingleCells') ,
		definition = function ( x, condition, files=F, A=NULL, B=NULL ) {
			if (!requireNamespace("MAST", quietly = TRUE)) {
				stop("MAST needed for this function to work. Please install it.",
						call. = FALSE)
			}
			toM <- function (x) {
				d <- as.matrix(x$dat)
				d[which(d==-20)] <- NA
				d[is.na(d)] <- 0
				d
			}
			if ( is.null(x$samples[,condition]) ) {
				stop( paste("the condition",condition, "is not defined in the samples table!"))
			}
			if ( is.null(A) ) {
				name = paste ("SingleCellAssay",condition)
				a <- x
			}else {
				keep <- which( x$samples[,condition] ==A | x$samples[,condition] == B)
				name = paste ("SingleCellAssay",condition, A, B)
				a <- reduceTo( x, what='col',to= colnames(x$data())[keep], name=name)
			}
			d <- toM(a)
			sca <- MAST::FromMatrix(class='SingleCellAssay', 
					exprsArray=d, 
					cData=data.frame(wellKey=colnames(d), GroupName = a$samples[,condition]), 
					fData=data.frame(primerid=rownames(d)))
			
			#groups <- sca@elementMetadata$GroupName <- a$samples[,condition]
			zlm.output <- MAST::zlm.SingleCellAssay(~ GroupName, sca, method='glm', ebayes=T)
			zlm.lr <- MAST::lrTest(zlm.output,'GroupName')
			
			x <- add_to_stat ( x, zlm.lr[,,'Pr(>Chisq)'], name )

			rm(sca)
			rm(zlm.output)
			rm(zlm.lr)
			gc(FALSE)
			invisible(x)
			
		}
)


add_to_stat <- function( x, stat, name ) {
	if ( length( match('p.adj BF',colnames(stat) )) == 0 & length( match('hurdle',colnames(stat) ))) {
		stat = cbind( stat, 'p.adj BF' = p.adjust(stat[,'hurdle'], method='bonferroni') )
	}
	if ( ! is.na( match( name, names(x$stats)))){
		x$stats[[ match( name, names(x$stats)) ]] <- stat
	}else {
		x$stats[[ length( x$stats ) +1 ]] <- stat
		names(x$stats)[length(x$stats) ] <- name
	}
	invisible(x)
}
