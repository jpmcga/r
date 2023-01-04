library(Seurat)
library(dplyr)
library(sctransform)
library(plyr)
library(magrittr)
library(tidyr)
library(openxlsx)


GemToSeuratObj <- function(gem.filepath = NULL,
                           min.tixels = 10,
                           mt.pattern = NULL){
  
  # GemToSeuratObj <- function(gem.filepath = NULL,
  #                            min.tixels = 10,
  #                            mt.pattern = "^mt."){ 
  # "^mt." for mouse, "^MT-" for human
  
  
  # Read in Gene Expression Matrix (GEM)
  gem <- read.table(gem.filepath,
                    sep = "\t",
                    header = TRUE,
                    row.names = 1)
  
  # Transpose GEM (columns = tixels, rows = genes)
  gem.transposed <- t(gem)
  
  # Create Seurat object
  seuratobj <- CreateSeuratObject(gem.transposed, min.cells = min.tixels)
  
  # Make percent.mt metadata column for normalizing with mitochondrial genes
  seuratobj <- PercentageFeatureSet(seuratobj,
                                    pattern =  mt.pattern, 
                                    col.name = "percent.mt") 
  
  # Normalize the data 
  seuratobj <- SCTransform(seuratobj,
                           vars.to.regress = "percent.mt", 
                           verbose = FALSE) 
  
  
  # Dimensionality reduction
  seuratobj <- RunPCA(seuratobj, verbose = FALSE)
  seuratobj <- RunUMAP(seuratobj, dims = 1:10, verbose = FALSE) 
  
  # Clustering
  seuratobj <- FindNeighbors(seuratobj, dims = 1:10, verbose = FALSE) 
  seuratobj <- FindClusters(seuratobj, resolution = 0.8, verbose = FALSE) 
  
  return(seuratobj)
}


GetTixelClusterDf <- function(seuratobj){
  
  tixel.idents <- Idents(seuratobj)
  tixel.idents.df <- data.frame(tixel.idents)
  tixel.idents.df <- data.frame(X = row.names(tixel.idents.df),
                                count = tixel.idents.df$tixel.idents)
  tixel.idents.df <- tixel.idents.df %>% separate(X,
                                                  c("x_cord", "y_cord"),
                                                  sep = "x")
  
  return(tixel.idents.df)
}


GetTixelClusterDf <- function(seuratobj){
  
  # Get factors (clusters) for each tixel and convert to df with columns for 
  # x and y coordinates.
  tixel.idents <- Idents(seuratobj)
  tixel.idents.df <- data.frame(tixel.idents)
  tixel.idents.df <- data.frame(X = row.names(tixel.idents.df),
                                count = tixel.idents.df$tixel.idents)
  tixel.idents.df <- tixel.idents.df %>% separate(X,
                                                  c("x_cord", "y_cord"),
                                                  sep = "x")
  
  return(tixel.idents.df)
}


GetDEGs <- function(seuratobj){
  
  # DEGs are identified between different clusters, which identifies
  # cluster-specific marker genes that facilitate labeling different cell
  # populations.  This identifies positive and negative markers in a single
  # cluster compared to all other cells-- but we're showing only positives
  
  # Finding Deferentially Expressed Feature markers to define clusters
  DEG.markers <- FindAllMarkers(seuratobj,
                                only.pos = TRUE,
                                min.pct = 0,
                                logfc.threshold = 0.01) 
  
  return(DEG.markers)
}


GetTopnDEGMarkers <- function(degs, n = 10){
  
  # Return top n markers for each cluster from full markers list
  
  markers <- degs %>%
    group_by(cluster) %>%
    top_n(n = n,
          wt = avg_log2FC) 
  return(markers)
}





