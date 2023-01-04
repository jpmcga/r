library(ggplot2)
library(dplyr)
library(plyr)
library(magrittr)
library(tidyr)
library(OpenImageR)
library(ggpubr)
library(grid)


PlotUMAPClusters <- function(seuratobj,
                             title = "UMAP Clustering"){
  
  DimPlot(seuratobj, label = TRUE,) +
    ggtitle(title) +  
    theme(plot.title = element_text(hjust = 0.5)) 
}


PlotSpatialClusters <- function(cluster.df,
                                tixel.size = 3,
                                title = "UMAP Spatial Clustering",
                                image.file = NULL){
  
  if (! is.null(image.file)){
    image <- OpenImageR::readImage(image.file)
    image.raster <- rasterGrob(image,
                               width = unit(1,"npc"),
                               height = unit(1,"npc"),
                               interpolate = FALSE)
  }
  
  ggplot(cluster.df, aes(x = as.numeric(x_cord),
                         y = as.numeric(y_cord),
                         color = count)) +
    ggtitle(title) + 
    {if (! is.null(image.file))
      annotation_custom(image.raster,
                        xmin=-Inf,
                        xmax=Inf,
                        ymin=-Inf,
                        ymax=Inf)
    } +
    geom_point(shape = 15,
               size = tixel.size) + 
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X",
                       limits = c(NA, NA), 
                       expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y",
                    limits = c(NA, NA), 
                    expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),
                ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5,
                                    size = 25,
                                    face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,
                                  face="bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}


PlotDEGTopnHeatMap <- function(seuratobj, deg.markers){
  
  heatmap <- DoHeatmap(seuratobj,
                       features = deg.markers$gene) +
    scale_fill_gradientn(colors = c("red", "black","green")) +
    theme(text = element_text(size = 6))
  
  return(heatmap)
}



