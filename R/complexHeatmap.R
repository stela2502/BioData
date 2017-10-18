
#' @name complexHeatmap
#' @aliases complexHeatmap,BioData-method
#' @rdname complexHeatmap-methods
#' @docType methods
#' @description plot the PCR heatmap using the heatmap.3 function included in this package 
#' @param x the BioData object
#' @param ofile the outfile to create in the x@outpath folder
#' @param colGroups columns in the samples table to use to order the data (first == order)
#' @param rowGroups rows in the annotation table to use to color the heatmap rows (first == order)
#' @param colColors a named list of column color vectors
#' @param rowColors a named list of row color vectors
#' @param pdf export as pdf (default = FALSE)
#' @param subpath the subpath for the plots (default = '')
#' @param heapmapCols the color function to calculate the heatmap colours ( default function (x) { c("darkgrey",bluered(x)) } )
#' @param brks how many breaks should the expression value color key have (default=10)
#' @title description of function complexHeatmap
#' @export 
