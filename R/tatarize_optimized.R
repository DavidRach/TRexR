#' This does this
#' 
#' @param num_colors Something
#' @param set_seed Default 16
#' 
#' @return Something
#' 
#' @export
#' 
#' @examples A <- 2 + 2
tatarize_optimized <- function(num_colors, set_seed=16) {
  set.seed(set_seed)
  colors_list = c()
  if (num_colors <= 24) {
    colors_list <- sample(tatarize_preferred, num_colors)
  } else if (num_colors <= 97) {
    colors_list <- sample(tatarize_preferred, 24)
    colors_list <- append(colors_list, sample(tatarize_secondary, (num_colors - 24))) 
  } else {
    colors_list <- sample(tatarize_preferred, 24)
    colors_list <- append(colors_list, sample(tatarize_secondary, 73))
    colors_list <- append(colors_list, sample(tatarize_avoid, (num_colors - 97)))
  }
  return(colors_list)
}