#' This does this
#' 
#' @param binned.data Something
#' @param export Default FALSE
#' @param outpath Default NULL
#' 
#' @importFrom utils write.csv
#' 
#' @return Something
#' 
#' @export
#' 
#' @examples A <- 2 + 2
TREX_results <- function(binned.data, 
                         export = FALSE,
                         outpath = NULL) {
  
  sample.table <- data.frame(total_cells = nrow(binned.data))
  sums <- vector()
  num_bins <- length(levels(binned.data$cuts))
  for (i in levels(binned.data$cuts)) {
    sums <- append(sums, sum(binned.data$cuts == i))
  }
  sample.table[, c(2:(num_bins + 1))] <- sums
  colnames(sample.table)[2:(num_bins + 1)] <- levels(binned.data$cuts)
  percent <- 100*(sums/nrow(binned.data))
  sample.table$degree_of_change <- (sum(percent[1] + percent[length(percent)]))
  sample.table$direction_of_change <- (sums[length(sums)] - sums[1]) / (sums[length(sums)] + sums[1])
  
  if (export) {

    if (is.null(outpath)){outpath <- getwd()}

    theFileName <- paste0(strftime(Sys.time(),"%Y-%m-%d_%H_%M"), "_TREX_results.csv")
    StashHere <- file.path(outpath, theFileName)

    write.csv(sample.table, StashHere, row.names = FALSE)
  }
  
  return(sample.table) 
}