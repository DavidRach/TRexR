# T-REX Analysis on AML Data

## Tracking Responders Expanding

``` r

knitr::opts_chunk$set(
  fig.width = 8,
  fig.height = 8,
  dpi = 96
)
```

``` r

Folder <- system.file("extdata", package="TRexR")
StorageLocation <- file.path(Folder, "aml")
fcs_files <- list.files(StorageLocation, pattern=".fcs", full.names=TRUE)

OutputLocation <- file.path(tempdir(), "TRexR-vignette")
dir.create(OutputLocation, showWarnings = FALSE, recursive = TRUE)
```

``` r

# read in files 
first.dataset <- as.data.frame(lapply(lapply(fcs_files[[1]], read.FCS), exprs), check.names=FALSE)
second.dataset <- as.data.frame(lapply(lapply(fcs_files[[2]], read.FCS), exprs), check.names=FALSE)

# save original column names for final export to FCS files 
orig.names <- colnames(first.dataset)

# rename columns with marker names 
colnames(first.dataset)[1:(length(first.dataset))] <- as.character(read.FCS(fcs_files[[1]])@parameters@data[["desc"]])
colnames(first.dataset)[1:2] <- c("Time", "Cell_length")
colnames(second.dataset)[1:(length(second.dataset))] <- as.character(read.FCS(fcs_files[[2]])@parameters@data[["desc"]])

# create meaningful names for data sets to be analyzed  
first.dataset$file_ID = "Timepoint 1"
second.dataset$file_ID = "Timepoint 2"
```

``` r

# find number of data points in smaller sample 
smaller_sample = nrow(first.dataset) 
if (nrow(second.dataset) < smaller_sample) { 
  smaller_sample <- nrow(second.dataset) 
}

# set a seed to get a reproducible sample 
set.seed(1)

# get an equal number of data points from both data sets
sampled.data = as.data.frame(rbind(
  first.dataset[sample(nrow(first.dataset), smaller_sample), ], 
  second.dataset[sample(nrow(second.dataset), smaller_sample), ]
))
```

``` r

colnames(sampled.data)
```

    ##  [1] "Time"        "Cell_length" "Barium"      "Barium"      "CD235a-141" 
    ##  [6] "CD19-142"    "CD117-143"   "CD11b-144"   "CD4-145"     "CD64-146"   
    ## [11] "CD7-147"     "CD34-148"    "CD61-150"    "CD123-151"   "CD13-152"   
    ## [16] "CD62L-153"   "CD45-154"    "CD183-156"   "CD33-158"    "CD11c-159"  
    ## [21] "CD14-160"    "CD15-164"    "CD16-165"    "CD24-166"    "CD38-167"   
    ## [26] "CD25-169"    "CD3-170"     "CD185-171"   "HLA-DR-174"  "CD184-175"  
    ## [31] "CD56-176"    "NA191"       "NA193"       "Lead"        "Lead"       
    ## [36] "Lead"        "file_ID"

``` r

# choose markers to make low dimensional embedding 
chosen.data = sampled.data[, c(5:31)]

# scale data 
transformed.data = as.data.frame(t(apply(chosen.data, 1, function(x) asinh(x/15))))
```

``` r

dataset.names = unique(sampled.data$file_ID)
#output.folder = paste0(getwd(), "/data_files/aml/T-REX output, ", dataset.names[1], " vs ", dataset.names[2])
#dir.create(output.folder)
```

``` r

# create UMAP (low dimensional embedding)
myumap <- umap(transformed.data, ret_model = TRUE, verbose = TRUE)
```

    ## 20:35:36 UMAP embedding parameters a = 1.896 b = 0.8006

    ## 20:35:36 Converting dataframe to numerical matrix

    ## 20:35:36 Read 63138 rows and found 27 numeric columns

    ## 20:35:36 Using Annoy for neighbor search, n_neighbors = 15

    ## 20:35:36 Building Annoy index with metric = euclidean, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 20:35:40 Writing NN index file to temp file /tmp/RtmpZwdUul/file54b67e8baaf7
    ## 20:35:40 Searching Annoy index using 2 threads, search_k = 1500
    ## 20:35:46 Annoy recall = 100%
    ## 20:35:46 Commencing smooth kNN distance calibration using 2 threads with target n_neighbors = 15
    ## 20:35:47 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 20:35:49 Commencing optimization for 200 epochs, with 1419788 positive edges
    ## 20:35:49 Using rng type: pcg
    ## 20:36:09 Optimization finished

``` r

# uncomment the line below to save the UMAP 
#save_uwot(myumap, paste0(strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),"_umap_model_results"))

umap.data = as.data.frame(myumap$embedding)
colnames(umap.data) <- c("UMAP1", "UMAP2")
```

``` r

my.binned = TREX(cbind(umap.data, file_ID = sampled.data$file_ID))

TREX_plot(
  binned.data = my.binned, 
  embed.type = "UMAP",
  caption = "Data from Ferrell et al., PLOS One. 2016",
  export = TRUE,
  outpath = OutputLocation
)
```

![](AML_files/figure-html/run%20TREX-1.png)

``` r

# The plot title will be incorrectly placed on the image below, but correctly placed on the exported png. The location is based on aspect ratio of the image, and can be adjusted with a title.height parameter if needed. 
```

``` r

TREX_results(my.binned, export = TRUE, outpath = OutputLocation)
```

    ##   total_cells [0,5] (5,15] (15,85) [85,95) [95,100] degree_of_change
    ## 1       63138 31499      9     130     107    31393         99.61038
    ##   direction_of_change
    ## 1        -0.001685429

``` r

my.clusters = TREX_cluster(
  binned.data = my.binned,
  db.eps = 0.3,
  marker.data = chosen.data
)

TREX_cluster_plot(
  cluster.data = my.clusters,
  binned.data = my.binned,
  embed.type = "UMAP",
  export = TRUE,
  outpath = OutputLocation
)
```

![](AML_files/figure-html/clustering-1.png)

``` r

TREX_cluster_results(my.clusters, export = TRUE, outpath = OutputLocation)
```

    ##    cluster median.change  mean.change
    ## 1        1      0.000000   0.05247879
    ## 2        2      0.000000   0.00000000
    ## 3        3      3.333333   2.75000000
    ## 4        4      0.000000   0.00000000
    ## 5        5    100.000000  99.95219358
    ## 6        6    100.000000  99.96576124
    ## 7        7    100.000000  99.82384824
    ## 8        8    100.000000 100.00000000
    ## 9        9    100.000000 100.00000000
    ## 10      10    100.000000 100.00000000
    ## 11      11    100.000000 100.00000000

``` r

# run MEM on DBSCAN clusters
MEM.input = my.clusters[, c(6:ncol(my.clusters))]
MEM.input <- MEM.input[order(MEM.input$cluster), ]
MEM.output = MEM(
  MEM.input, 
  transform = TRUE, 
  cofactor = 15,
  zero.ref = TRUE, 
  new.marker.names = "CD235a,CD19,CD117,CD11b,CD4,CD64,CD7,CD34,CD61,CD123,CD13,CD62L,CD45,CD183,CD33,CD11c,CD14,CD15,CD16,CD24,CD38,CD25,CD3,CD185,HLA-DR,CD184,CD56", 
  scale.matrix = "arcsinh",
  scale.factor = 1
)

build_heatmaps(
  MEM.output, 
  cluster.MEM = "none",
  output.files = TRUE, # No outpath arguments 
)
```

![](AML_files/figure-html/MEM-1.png)![](AML_files/figure-html/MEM-2.png)![](AML_files/figure-html/MEM-3.png)

    ##  [1] "MEM label for cluster 1 : ▲ CD64+4 CD45+4 CD33+4 HLA-DR+4 CD11b+3 CD4+3 CD184+3 CD11c+2 CD38+2 CD183+1 CD56+1 ▼ None"                                                                        
    ##  [2] "MEM label for cluster 2 : ▲ CD235a+8 CD19+3 CD33+3 CD4+1 CD45+1 CD183+1 HLA-DR+1 ▼ None"                                                                                                     
    ##  [3] "MEM label for cluster 3 : ▲ CD45+8 CD3+7 CD4+6 CD64+5 CD7+5 HLA-DR+5 CD184+5 CD183+4 CD33+4 CD38+4 CD11b+3 CD56+3 CD11c+2 CD185+2 CD34+1 CD123+1 CD13+1 CD62L+1 ▼ None"                      
    ##  [4] "MEM label for cluster 4 : ▲ CD45+8 CD7+7 CD16+5 CD56+5 CD11b+4 CD38+4 CD64+3 CD183+3 CD184+3 CD34+2 CD33+2 CD11c+2 HLA-DR+2 CD4+1 CD13+1 CD14+1 ▼ None"                                      
    ##  [5] "MEM label for cluster 5 : ▲ CD3+10 CD45+9 CD7+6 CD184+6 CD4+4 CD62L+4 CD34+2 CD56+2 CD64+1 CD13+1 CD183+1 ▼ None"                                                                            
    ##  [6] "MEM label for cluster 6 : ▲ CD3+10 CD45+8 HLA-DR+5 CD184+5 CD11b+4 CD4+2 CD64+2 CD33+2 CD11c+2 CD56+2 CD19+1 CD61+1 CD123+1 CD13+1 CD62L+1 CD14+1 CD38+1 CD185+1 ▼ None"                     
    ##  [7] "MEM label for cluster 7 : ▲ CD3+10 CD64+9 CD45+8 CD7+6 CD184+5 CD4+4 CD62L+3 CD56+2 CD235a+1 CD19+1 CD11b+1 CD34+1 HLA-DR+1 ▼ None"                                                          
    ##  [8] "MEM label for cluster 8 : ▲ CD3+10 CD235a+8 CD45+8 HLA-DR+6 CD19+5 CD64+5 CD184+5 CD62L+4 CD56+4 CD14+3 CD185+3 CD7+2 CD34+2 CD183+2 CD4+1 CD61+1 CD13+1 CD11c+1 CD15+1 CD24+1 CD38+1 ▼ None"
    ##  [9] "MEM label for cluster 9 : ▲ CD3+10 CD235a+8 CD45+8 CD7+6 CD184+6 CD4+4 CD62L+4 CD19+3 CD64+3 CD34+3 CD56+3 CD183+2 CD117+1 CD61+1 CD13+1 CD14+1 CD15+1 HLA-DR+1 ▼ None"                      
    ## [10] "MEM label for cluster 10 : ▲ CD3+10 CD34+9 CD45+7 CD19+6 CD64+6 CD7+6 CD4+5 CD15+5 CD235a+4 CD11b+4 CD184+4 CD117+3 CD62L+1 CD56+1 ▼ None"                                                   
    ## [11] "MEM label for cluster 11 : ▲ CD3+10 CD45+8 CD11b+7 CD62L+6 CD64+4 CD61+3 CD14+3 CD15+3 CD184+3 CD7+2 CD34+2 HLA-DR+2 CD19+1 CD13+1 CD183+1 CD11c+1 CD24+1 CD185+1 ▼ None"

``` r

TREX_counts(my.clusters, export = TRUE, outpath = OutputLocation)
```

    ##    cluster total_cells [0,5] [95,100]
    ## 1        1       31346 31346        0
    ## 2        2          65    65        0
    ## 3        3          60    60        0
    ## 4        4          28    28        0
    ## 5        5       17536     0    17536
    ## 6        6       13143     0    13143
    ## 7        7         615     0      615
    ## 8        8          14     0       14
    ## 9        9          71     0       71
    ## 10      10          13     0       13
    ## 11      11           1     0        1

``` r

# create heatmaps of every marker across UMAP axes 
umap.bychannel <- as_tibble(umap.data) %>%
  bind_cols(transformed.data)  %>%
  gather(channel, intensity, -UMAP1, -UMAP2) %>%
  mutate(across(channel,factor))%>%
  group_split(channel) %>%
  map(
    ~ggplot(., aes(x = UMAP1, y = UMAP2, col = intensity)) +
      geom_point(shape = ".", size = 8) +
      scale_color_gradientn(
        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(5)) +
      facet_grid(~ channel, labeller = function(x) label_value(x, multi_line = FALSE)) +
      coord_fixed() +
      theme_bw() +
      theme(
        strip.text.x = element_text(size = 20),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    ) %>%
  plot_grid(plotlist = ., align = 'hv', ncol = 8)

theFileName <- paste0(strftime(Sys.time(),"%Y-%m-%d_%H-%M"),"UMAP on transformed data.png")
StashHere <- file.path(OutputLocation, theFileName)

png(
  StashHere,
  height = 2000,
  width = 4000
)
print(umap.bychannel)
dev.off()
```

    ## agg_png 
    ##       2

``` r

# export new FCS files with sampled data, percent change, and umap axes
to.export = cbind(sampled.data, umap.data, percent.change = my.binned$percent.change)

# add column for cluster number 
to.export$cluster <- c(0)
bins.of.interest = unique(my.clusters$cuts)
temp.cuts = wafflecut(to.export$percent.change, levels(my.binned$cuts))
to.export[temp.cuts %in% bins.of.interest, ]$cluster <- my.clusters$cluster

desc = colnames(to.export)
desc <- desc[! desc %in% c("file_ID")]
colnames(to.export)[1:ncol(sampled.data) - 1] <- orig.names
sep.files = split(to.export, to.export$`file_ID`)
for (i in 1:length(sep.files)){
  reduce.data = subset(sep.files[[i]], select = -c(`file_ID`))
  mat.input <- as.matrix(reduce.data)
  metadata <- data.frame(name = dimnames(mat.input)[[2]], desc = desc)
  metadata$range <- apply(apply(mat.input, 2, range), 2, diff)
  metadata$minRange <- apply(mat.input, 2, min)
  metadata$maxRange <- apply(mat.input, 2, max)
  input.flowframe <- new("flowFrame", exprs = mat.input, parameters = AnnotatedDataFrame(metadata))
  newname = str_remove(basename(fcs_files[i]), ".fcs")
  new.filename = file.path(OutputLocation, paste0(strftime(Sys.time(),"%Y-%m-%d_%H-%M"), "_", newname, "_T-REX.fcs"))
  write.FCS(input.flowframe, filename = new.filename)
  print(paste("FCS file ", i, " done", sep = ""))
}
```

    ## [1] "FCS file 1 done"
    ## [1] "FCS file 2 done"

``` r

# MEM Clean Up, still making an output files folder, remove it. 
if (dir.exists("output files")){unlink("output files", recursive=TRUE)}
```

``` r

sessionInfo()
```

    ## R version 4.6.2 (2026-10-29)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] TRexR_0.99.1        fancycut_0.1.4      cowplot_1.2.0      
    ##  [4] purrr_1.2.2         stringr_1.6.0       dplyr_1.2.1        
    ##  [7] RColorBrewer_1.1-3  Biobase_2.72.0      BiocGenerics_0.58.1
    ## [10] generics_0.1.4      tidyr_1.3.2         cytoMEM_1.16.0     
    ## [13] ggplot2_4.0.3       uwot_0.2.5          Matrix_1.7-6       
    ## [16] flowCore_2.24.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6        xfun_0.62           bslib_0.12.0       
    ##  [4] caTools_1.18.4      lattice_0.23-1      vctrs_0.7.3        
    ##  [7] tools_4.6.2         bitops_1.1-0        stats4_4.6.2       
    ## [10] tibble_3.3.1        pkgconfig_2.0.3     KernSmooth_2.23-27 
    ## [13] S7_0.2.2            desc_1.4.3          S4Vectors_0.50.3   
    ## [16] lifecycle_1.0.5     compiler_4.6.2      farver_2.1.2       
    ## [19] FNN_1.1.4.1         textshaping_1.0.5   gplots_3.3.0       
    ## [22] codetools_0.2-20    htmltools_0.5.9     sass_0.4.10        
    ## [25] yaml_2.3.12         pillar_1.11.1       pkgdown_2.2.1      
    ## [28] jquerylib_0.1.4     cachem_1.1.0        dbscan_1.2.7       
    ## [31] RSpectra_0.16-2     RProtoBufLib_2.24.0 gtools_3.9.5       
    ## [34] tidyselect_1.2.1    digest_0.6.39       stringi_1.8.9      
    ## [37] labeling_0.4.3      fastmap_1.2.0       grid_4.6.2         
    ## [40] cli_3.6.6           magrittr_2.0.5      withr_3.0.3        
    ## [43] scales_1.4.0        rmarkdown_2.32      matrixStats_1.5.0  
    ## [46] otel_0.3.0          gridExtra_2.3.1     cytolib_2.24.0     
    ## [49] ragg_1.5.2          evaluate_1.0.5      knitr_1.52         
    ## [52] RcppAnnoy_0.0.23    rlang_1.3.0         Rcpp_1.1.2         
    ## [55] glue_1.8.1          jsonlite_2.0.0      R6_2.6.1           
    ## [58] systemfonts_1.3.2   fs_2.1.0
