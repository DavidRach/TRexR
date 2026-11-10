#' This does this
#' 
#' @param binned.data Something
#' @param bins.of.interest Default NULL
#' @param db.eps Default 0.5
#' @param marker.data Default NULL
#' 
#' @importFrom dplyr filter %>%
#'  
#' @return Something
#' 
#' @export
#' 
#' @examples A <- 2 + 2
TREX_cluster <- function(binned.data, 
                         bins.of.interest = NULL,
                         db.eps = 0.5,
                         marker.data = NULL) {           
  
  if (!is.null(marker.data)) {
    if ("file_ID" %in% colnames(marker.data)) {
      marker.data <- subset(marker.data, select = -c(file_ID))
    }
    binned.data <- cbind(binned.data, marker.data)
  }
  
  if (is.null(bins.of.interest)) {
    bin.levels = levels(binned.data$cuts)
    bins.of.interest <- c(bin.levels[1], bin.levels[length(bin.levels)])
  }  
  
  regions.of.interest = binned.data[binned.data$cuts %in% bins.of.interest, ]
  db.result = dbscan::dbscan(regions.of.interest[, 1:2], eps = db.eps, minPts = 1)
  cluster.data = cbind(regions.of.interest, cluster = db.result$cluster)
  cluster.data <- cluster.data %>%
    filter(cluster != 0)
  
  return(cluster.data)
}