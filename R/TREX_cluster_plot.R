#' This does this
#' 
#' @param cluster.data Something
#' @param binned.data Default NULL
#' @param embed.type Default "Embedding"
#' @param colors Default NULL
#' @param export Default FALSE
#' @param outpath Default NULL
#' 
#' @importFrom ggplot2 ggplot geom_point aes coord_fixed scale_color_manual
#'  labs guides guide_legend ggsave
#' 
#' @return Something
#' 
#' @export
#' 
#' @examples A <- 2 + 2
TREX_cluster_plot <- function(cluster.data,
                              binned.data = NULL,
                              embed.type = "Embedding",
                              colors = NULL,
                              export = FALSE,
                              outpath = NULL) {
  
  embed.x = paste(embed.type, "1")
  embed.y = paste(embed.type, "2")
  
  cluster.plot = ggplot()
  color.values = tatarize_optimized(length(unique(as.factor(cluster.data$cluster))))
  
  if (!is.null(colors)) {
    color.values <- colors
  } 
  
  # use all embedding points to calculate ratio if included 
  if (!is.null(binned.data)) {
    cluster.plot <- cluster.plot + geom_point(data = binned.data, aes(x = x, y = y), col = "lightgray")
    color.values <- append(color.values, "lightgray")
    range <- apply(apply(binned.data[, 1:2], 2, range), 2, diff)
    graphical.ratio <- (range[1] / range[2])
  } else {
    range <- apply(apply(cluster.data[, 1:2], 2, range), 2, diff)
    graphical.ratio <- (range[1] / range[2])
  }
  
  cluster.plot <- cluster.plot + 
    geom_point(data = cluster.data, aes(x = x, y = y, col = as.factor(cluster)), cex = 1.5) +
    coord_fixed(ratio = graphical.ratio) +
    scale_color_manual(values = color.values) +
    labs(
      x = embed.x, 
      y = embed.y, 
      title = "DBSCAN Clusters, bins of interest"
    ) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme_TREX()
  
  if (export) {

    if (is.null(outpath)){outpath <- getwd()}
    
    ggsave(
      filename = paste0(strftime(Sys.time(), "%Y-%m-%d_%H_%M"), "_cluster_plot.png"), 
      plot = cluster.plot,
      path     = outpath,
      width = 8, 
      height = 8
    )
  }
  
  return(cluster.plot)
}