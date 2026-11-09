#' This does this
#' 
#' @param bins Something
#' 
#' @importFrom grDevices colorRampPalette
#' 
#' @return Something
#' 
#' @export
get_TREX_colors <- function(bins) {
  if (length(bins) == 3) {
    bin.colors = default_colors[c(1,3,5)]
  } else if (length(bins) == 4) {
    bin.colors = default_colors[c(1,2,4,5)]
  } else if (length(bins) == 5) {
    bin.colors = default_colors
  } else if(length(bins) > 5 & length(bins) < 13) {
    bin.colors = colorRampPalette(default_colors)(length(bins))
  } else if(length(bins) == 13) {
    bin.colors = thirteen_colors
  } else if (length(bins) > 13) {
    # library(RColorBrewer)
    bin.colors = colorRampPalette(thirteen_colors)(length(bins))
  }
  
  names(bin.colors) <- bins
  return(bin.colors)
}