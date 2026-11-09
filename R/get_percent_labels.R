
#' This does this
#' 
#' @param bin.levels Something
#' @param data.names Something
#' 
#' @importFrom stringr str_detect str_extract
#' 
#' @return Something
#' 
#' @export
get_percent_labels <- function(bin.levels, data.names) {
  bin.labels = vector()
  for (bin in bin.levels) {
    num.list = strsplit(bin, split = ",")
    first.num = num.list[[1]][1]
    second.num = num.list[[1]][2]
    
    # check for decimal points 
    if (str_detect(first.num, "\\.")) {
      first.num <- as.numeric(str_extract(first.num, "[0-9]+\\.[0-9]+"))
    } else {
      first.num <- as.numeric(str_extract(first.num, "[0-9]+"))
    }
    if (str_detect(second.num, "\\.")) {
      second.num <- as.numeric(str_extract(second.num, "[0-9]+\\.[0-9]+"))
    } else {
      second.num <- as.numeric(str_extract(second.num, "[0-9]+"))
    }
    
    if (first.num == 0) {
      cluster.label = paste0("\u2265", (100 - second.num), "% from ", data.names[1])
    } else if (second.num == 100) {
      cluster.label = paste0("\u2265", first.num, "% from ", data.names[2])
    } else if (second.num <= 50 & first.num <= 50) {
      cluster.label = paste0((100 - second.num), "-", (100 - first.num), "% from ", data.names[1]) 
    } else if (second.num > 50 & first.num > 50) {
      cluster.label = paste0(first.num, "-", second.num, "% from ", data.names[2])
    } else {
      cluster.label = paste0("from ", data.names[1], " and ", data.names[2])
    }
    bin.labels <- append(bin.labels, cluster.label)
  }
  return(bin.labels)
}