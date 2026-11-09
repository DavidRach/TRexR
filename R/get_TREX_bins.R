#' This does this
#' 
#' @param bin.size Something
#' 
#' @return Something
#' 
#' @export
get_TREX_bins <- function(bin.size) {
  if (bin.size >= 100 | bin.size <= 0) {
    stop("To generate bins from a percentage or number, the value must be between 1 and 100")
  }
  
  n_bins = round(100/bin.size, digits = 0)
  bin.vector = vector("character", n_bins)
  start.b = "("
  end.b = "]"
  bin.vector[1] <- paste0("[0,",bin.size,"]")
  center.bin = floor(n_bins/2)
  
  # handle special case of only 3 bins 
  if (center.bin == 1) {
    center.bin <- 2 
  }
  
  for(i in 2:n_bins) {
    lower.num = (i - 1)*bin.size
    upper.num = i*bin.size 
    
    # last bin in list
    if(i == n_bins) {              
      bin.vector[i] <- paste0("[",lower.num,",100]") 
      
      # center bin
    } else if(i == center.bin) {   
      bin.vector[i] <- paste0("(",lower.num,",",upper.num,")") 
      start.b <- "["
      end.b <- ")"
      
      # all other bins
    } else {                            
      bin.vector[i] <- paste0(start.b,lower.num,",",upper.num,end.b)         
    }
  } 
  
  return(bin.vector)
}