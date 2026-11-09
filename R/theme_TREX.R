#' This does this
#' 
#' @importFrom ggplot2 theme_bw theme element_blank element_text %+replace%
#' 
#' @return Something
#' 
#' @export
theme_TREX <- function() {
  theme_bw() %+replace% 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 14),
      legend.title = element_blank()
    )
}