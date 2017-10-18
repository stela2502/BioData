#' @name extractCodonInformation
#' @aliases extractCodonInformation,tRNAMINT-method
#' @rdname extractCodonInformation-methods
#' @docType methods
#' @description Converts the MINT tRNA annoation strings into codon information stored in the annotation.
#' @param x the tRNAMINT object
#' @param col the column conatining the codon information default= "Sequence.locations.in.tRNA.space..comma.deliminated."
#' @title description of function extractCodonInformation
#' @export 
if ( ! isGeneric('extractCodonInformation') ){ setGeneric('extractCodonInformation', ## Name
		function ( x, col= "Sequence.locations.in.tRNA.space..comma.deliminated." ) { ## Argumente der generischen Funktion
			standardGeneric('extractCodonInformation') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)
}else {
	print ("Onload warn generic function 'extractCodonInformation' already defined - no overloading here!")
}

setMethod('extractCodonInformation', signature = c ('tRNAMINT'),
		definition = function ( x, col= "Sequence.locations.in.tRNA.space..comma.deliminated." ) {
			# get a vector of Sequence.locations.in.tRNA.space..comma.deliminated.
			# and returns a list of Sequence.locations.in.tRNA.space per vector entry
			split_to_tRNAs <- function( x) {
				stringr::str_split( x ,', ')
			}
			# getsone list entry of the list returned from split_to_tRNAs
			# which is a vector of Sequence.locations.in.tRNA.space
			# it then has to ectract the codon information
			# trna133_GlyCCC_1_-_16872434_16872504@1.52.52
			split_to_codon <- function(x) {
				sort(unique(unlist(lapply( stringr::str_split(x,'_'), function( a) { a[2]} ))))
			}
			codon_list <- lapply( split_to_tRNAs(  x$annotation[,col] ), split_to_codon )
			
			x$usedObj$Codons <- names(table(unlist(codon_list)))
			
			for ( codon in names(table(unlist(codon_list)))) {
				x$annotation[,codon] = 0
			}
			for ( i in 1:length(codon_list)){
				for (codon in codon_list[[i]] ) {
					x$annotation[i,codon] = 1
				}
			}
			# now all codons are stored in the annotation table
			# start with the sample summaries
			invisible(x)
		}
)