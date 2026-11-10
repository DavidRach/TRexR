#' This does this
#' 
#' @param cluster.data Something
#' @param export Default FALSE
#' @param outpath Default NULL
#' 
#' @importFrom stats median
#' @importFrom utils write.csv
#' 
#' @return Something
#' 
#' @export
#' 
#' @examples A <- 2 + 2
TREX_cluster_results <- function(cluster.data,
                                 export = FALSE,
                                 outpath = NULL) {
  
  split.data = split(cluster.data, cluster.data$cluster)
  median.percent.change = lapply(split.data, function(x) median(x[, which(colnames(x) == "percent.change")]))
  mean.percent.change = lapply(split.data, function(x) mean(x[, which(colnames(x) == "percent.change")]))
  
  results = data.frame(
    cluster = unique(cluster.data$cluster),
    median.change = unlist(median.percent.change),
    mean.change = unlist(mean.percent.change)
  )
  
  if (export) {
    if (is.null(outpath)){outpath <- getwd()}
    theFileName <- paste0(strftime(Sys.time(),"%Y-%m-%d_%H%M"),"_cluster_percent_change.csv")
    StashHere <- file.path(outpath, theFileName)
    write.csv(mean.percent.change, StashHere, row.names = FALSE)
  }
  
  return(results)
}