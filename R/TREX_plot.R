#' This does this
#' 
#' @param binned.data Something
#' @param title.height Default -3
#' @param embed.type Default "Embedding"
#' @param percent.labels Default TRUE
#' @param caption Default NULL
#' @param return.obj Default "drawn plot"
#' @param export Default FALSE
#' @param outpath Default NULL
#' 
#' @importFrom ggplot2 ggplot geom_point aes coord_fixed scale_color_manual
#'  guides guide_legend labs waiver ggsave
#' @importFrom grid grobTree gpar textGrob unit grobWidth
#' @importFrom gridExtra arrangeGrob
#' @importFrom cowplot ggdraw
#' 
#' @return Something
#' 
#' @export
TREX_plot <- function(binned.data,
                      title.height = -3,
                      embed.type = "Embedding",
                      percent.labels = TRUE,
                      caption = NULL,
                      return.obj = "drawn plot",
                      export = FALSE,
                      outpath = NULL) {

  # get dataset names from file_ID column 
  data.names = unique(binned.data$file_ID)
  
  # reorder bins so that regions of interest are on top 
  bins = levels(binned.data$cuts)
  if (length(bins) %% 2 == 0) {
    center.indx = length(bins)/2
    bin.levels = c(bins[center.indx], bins[center.indx + 1]) 
    for (i in 1:(center.indx - 1)) {
      bin.levels <- append(bin.levels, bins[center.indx - i])
      bin.levels <- append(bin.levels, bins[center.indx + i + 1])
    }
  } else {
    center.indx = round(length(bins)/2, 0) + 1
    bin.levels = bins[center.indx]
    for (i in 1:(center.indx - 1)) {
      bin.levels <- append(bin.levels, bins[center.indx - i])
      bin.levels <- append(bin.levels, bins[center.indx + i])
    }
  }
  
  binned.data$cuts <- factor(binned.data$cuts, levels = bin.levels)
  binned.data <- binned.data[order(binned.data$cuts), ]
  range <- apply(apply(binned.data[, 1:2], 2, range), 2, diff)
  graphical.ratio <- (range[1] / range[2])
  
  if (percent.labels) {
    bin.labels = get_percent_labels(bins, data.names)
  } else {
    bin.labels = waiver()
  }
  
  embed.x = paste(embed.type, "1")
  embed.y = paste(embed.type, "2")
  
  trex.plot = ggplot(binned.data) + 
    geom_point(aes(x = x, y = y, color = cuts), cex = 1) + 
    coord_fixed(ratio = graphical.ratio) +
    scale_color_manual(
      labels = bin.labels,
      values = get_TREX_colors(bins),
      breaks = bins
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(x =  embed.x, y = embed.y, caption = caption) +
    theme_TREX()
  
  # add red/blue colored title 
  titleGrobs <- grobTree(
    gp = gpar(fontsize = 16, fontface = "bold"),
    textGrob(label = data.names[1], 
             name = "title1",
             x = unit(3, "lines"), 
             y = unit(title.height, "lines"), 
             hjust = 0, vjust = 0, 
             gp = gpar(col = "#000080")),
    textGrob(label = " vs. ", 
             name = "title2",
             x = grobWidth("title1") + unit(3, "lines"), 
             y = unit(title.height, "lines"),
             hjust = 0, vjust = 0),
    textGrob(label = data.names[2], 
             name = "title3",
             x = grobWidth("title1") + grobWidth("title2") + unit(3, "lines"), 
             y = unit(title.height, "lines"),
             hjust = 0, vjust = 0, 
             gp = gpar(col = "#8B0000"))
  )
  trex.titled <- arrangeGrob(trex.plot, top = titleGrobs, padding = unit(2.6, "line"))
  
  if (export) {

    if (is.null(outpath)){outpath <- getwd()}

    ggsave(
      paste0(
        strftime(Sys.time(), "%Y-%m-%d_%H_%M "),
        data.names[1], " vs ", data.names[2],
        " T-REX plot.png"
      ), 
      plot = trex.titled,
      path = outpath, 
      width = 8, 
      height = 8
    )
  }
  
  if (return.obj == "drawn plot") { 
    return(ggdraw(trex.titled))
  } else if (return.obj == "ggplot") {
    return(trex.plot)
  } else if (return.obj == "gtable") {
    return(trex.titled)
  } else {
    stop("uninterpretable input for return.obj")
  }
  
}