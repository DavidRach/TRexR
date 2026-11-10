#' This does this
#' 
#' @param embedding.data Something
#' @param kvalue Default 60
#' @param bins Default c('[0,5]','(5,15]','(15,85)','[85,95)','[95,100]'))
#' 
#' @importFrom FNN knnx.index
#' @importFrom fancycut wafflecut
#' 
#' @return Something
#' 
#' @export
TREX <- function(embedding.data,                   
                 kvalue = 60,                       
                 bins = c('[0,5]','(5,15]','(15,85)','[85,95)','[95,100]')           
                 ) {
  
  # 'bins' here refers to groups created by T-REX based on % change. 
  # This variable controls the red and blue areas on the plot. 
  
  # warnings for improperly formatted data 
  if (is.na(match("file_ID", colnames(embedding.data)))) {
    stop("Embedding.data has no column named file_ID.")
  }

  data.name = unique(embedding.data$file_ID)[1]
  if (length(which(embedding.data$file_ID == data.name)) != nrow(embedding.data)/2) {
    stop("Datasets are not equally sampled according to file_ID.")
  }
  
  # create bins vector if given integer or % 
  if (is.numeric(bins)) {
    bins <- get_TREX_bins(round(100/bins, digits = 1))
  } else if(is.character(bins) & length(bins) == 1) {
    if(grepl("%", bins)) { 
      bin.size = as.numeric(substr(bins, start = 1, stop = nchar(bins) - 1))
      if (bin.size > 30) {
        stop("Cannot have less than three bins.")
      }
      bins <- get_TREX_bins(bin.size)
    }
  }
  
  # KNN search per cell 
  neighbor.index = knnx.index(embedding.data[, 1:2], embedding.data[, 1:2], k = kvalue)
  
  # assign dataset belonging by row number
  neighbor.index[neighbor.index <= nrow(embedding.data)/2] <- 0
  neighbor.index[neighbor.index > nrow(embedding.data)/2] <- 1

  # calculate percent change in each KNN region
  percent.change = (rowSums(neighbor.index) / kvalue * 100)
  
  # binning 
  binned.data <- data.frame(
    x = embedding.data[, 1],
    y = embedding.data[, 2],
    file_ID = embedding.data$file_ID,
    percent.change = percent.change
  )
  binned.data$cuts <- wafflecut(round(binned.data$percent.change), bins)
  
  return(binned.data)
}