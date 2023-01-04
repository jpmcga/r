library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(ggpubr)
library(grid)
library(cowplot)
require(scales)


# UMAP clustering with sample and reference as separate clusters
PlotSeparateClusters <- function(sample.ref.clustered){
  
  DimPlot(sample.ref.clustered, 
          reduction = "umap",
          group.by = "tech", 
          order = c("gem.seurat.obj","ref.seurat.obj") ,
          pt.size = 0.5,
          cols = c("blue", "magenta"))+
    ggtitle("Clustering, Sample vs Reference data")
}

# UMAP clustering with sample and reference data as combined clusters
PlotCombinedClusters <- function(sample.ref.clustered){
  
  DimPlot(sample.ref.clustered, 
          reduction = "umap", 
          label = TRUE,
          pt.size = 0.5)+
    ggtitle("Clustering, Sample and Reference data integrated")
  
}

# Cell annotation for the integration of sample and reference data with cell types
PlotCombinedCellAnn <- function(sample.ref.clustered){
  
  sample.ref.clustered$annotation = cells.under.10
  
  DimPlot(sample.ref.clustered, 
          reduction = "umap",
          group.by = "annotation", 
          pt.size = 1.2, 
          label = TRUE) + 
    ggtitle("Cell Annotation, integrated data")
  
}

# Cell annotation for the integration of sample and reference data with cell 
# types but with no cell types removed
PlotCombinedCellAnn2 <- function(sample.ref.clustered){
  
  sample.ref.clustered$annotation = cell.annotation
  
  DimPlot(sample.ref.clustered, 
          reduction = "umap",
          group.by = "annotation", 
          pt.size = 1.2, 
          label = TRUE) + 
    NoLegend() + 
    ggtitle("Cell Annotation, integrated data")
}

PlotSideBySide <- function(plot1, plot2){
  
  side.by.side.plots <- plot_grid(plot1,plot2)
  return(side.by.side.plots)
  
}



