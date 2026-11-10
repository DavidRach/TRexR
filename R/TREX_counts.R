#' This does this
#' 
#' @param cluster.data Something
#' @param export Default TRUE
#' @param outpath Default NULL
#' 
#' @importFrom utils write.csv
#' 
#' @return Something
#' 
#' @export
TREX_counts <- function(cluster.data,
                        export = TRUE,
                        outpath = NULL) {
  
  results = data.frame(cluster = unique(cluster.data$cluster))
  
  # count total cells in each cluster
  split.clusters = split(cluster.data, as.factor(cluster.data$cluster))
  results$total_cells <- sapply(split.clusters, NROW)
  
  # find bins of interest
  bins.of.interest = vector()
  for(i in split.clusters) {
    for (z in unique(i$cuts)) {
      if (!z %in% bins.of.interest) {
        bins.of.interest <- append(bins.of.interest, z)
      }
    }
  }
  bins.of.interest <- bins.of.interest[order(match(bins.of.interest, as.vector(levels(cluster.data$cuts))))]
  
  # count how many cells are in each bin 
  results <- cbind(results, data.frame = (matrix(0, ncol = length(bins.of.interest))))
  colnames(results)[3:ncol(results)] <- bins.of.interest
  
  col_indx = 3
  for (w in bins.of.interest) {
    for (i in 1:length(split.clusters)) {
      results[i, col_indx] <- sum(split.clusters[[i]][["cuts"]] == w)
    }
    col_indx <- col_indx + 1
  }
  
  if (export) {

    if (is.null(outpath)){outpath <- getwd()}
    theFileName <- paste(strftime(Sys.time(),"%Y-%m-%d_%H%M"),"_cluster_counts.csv")
    StashHere <- file.path(outpath, theFileName)
    write.csv(results, StashHere, row.names = FALSE)
  }
  
  return(results)
}