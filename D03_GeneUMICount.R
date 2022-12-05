library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

data <- read.table("D03_stdata_names.tsv", header=T, row.names = 1, as.is = TRUE)

object <- CreateSeuratObject(counts = data, assay = assay)

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), 
                       filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

brain <- object

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.5) + NoLegend()
plot2 <- VlnPlot(brain, features = "nFeature_Spatial", pt.size = 0.5) + NoLegend()
wrap_plots(plot1, plot2)

plot3 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
plot3$layers[[1]]$aes_params <- c(plot3$layers[[1]]$aes_params, shape=22) 
plot4 <- SpatialFeaturePlot(brain, features = "nFeature_Spatial") + theme(legend.position = "right")
plot4$layers[[1]]$aes_params <- c(plot4$layers[[1]]$aes_params, shape=22) 
wrap_plots(plot3, plot4)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
ElbowPlot(brain, ndims = 50)

brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE, resolution = 0.2)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p1
p2 <- SpatialDimPlot(brain, label = FALSE, label.size = 3, pt.size.factor = 2.5)
p2$layers[[1]]$aes_params <- c(p2$layers[[1]]$aes_params, shape=22) 
p2

wrap_plots(p1,p2)

de_markers <- FindMarkers(brain, ident.1 = 3)
p3 <- SpatialFeaturePlot(object = brain, features = rownames(de_markers)[2], alpha = c(0.1, 1), ncol = 1)
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
p3

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000], selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
p4 <- SpatialFeaturePlot(brain, features = top.features[6], ncol = 1, alpha = c(0.1, 1))
p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
p4

