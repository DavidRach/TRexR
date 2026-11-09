#' This does this
#' 
#' @param cluster.data Something
#' @param export Default FALSE
#' 
#' @importFrom stats median
#' @importFrom utils write.csv
#' 
#' @return Something
#' 
#' @export
TREX_cluster_results <- function(cluster.data,
                                 export = FALSE) {
  
  split.data = split(cluster.data, cluster.data$cluster)
  median.percent.change = lapply(split.data, function(x) median(x[, which(colnames(x) == "percent.change")]))
  mean.percent.change = lapply(split.data, function(x) mean(x[, which(colnames(x) == "percent.change")]))
  
  results = data.frame(
    cluster = unique(cluster.data$cluster),
    median.change = unlist(median.percent.change),
    mean.change = unlist(mean.percent.change)
  )
  
  if (export) {
    write.csv(mean.percent.change, paste0(strftime(Sys.time(),"%Y-%m-%d_%H%M"),"_cluster_percent_change.csv"), row.names = FALSE)
  }
  
  return(results)
}