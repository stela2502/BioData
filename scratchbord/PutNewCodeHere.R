
log <- function (x) {
	if ( ! x$logged ) {
		x$data <- apply( x$data, 2, function(x){ log(x+1) } )
	}
	invisible(x)
}