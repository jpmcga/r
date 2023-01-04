library(Seurat)
library(patchwork)
library(Matrix)
library(sctransform)
library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(SingleR)
library("monocle")
library("BiocManager")
library("SummarizedExperiment")
library("scater")
library("S4Vectors")
require(scales)


ProcessGEM <- function(gem.filepath){
  # Create Seurat object
  seurat.obj <- CreateSeuratObject(gem.filepath, min.cells = 10)
  
  return(seurat.obj)
}

FindIDNames <- function(ref.data){
  
  fData(ref.data) #featureData
  pData(ref.data) #access phenotypic data
  # use these to retrieve the expression matrix information (because the format 
  # is different from the format of the GEM which is in columns and rows)
  
  ref.df = data.frame(gene = fData(ref.data)$gene_short_name)
  idname = names(table(pData(ref.data)$id))
  embryo.ages = data.frame(pData(ref.data)$id, pData(ref.data)$day)
  embryo.ages <- unique(embryo.ages)
  # builds a table that we can use to identify which day of the embryo we're 
  # looking at by showing which day represents which ID)
  
  return(embryo.ages)
}


BuildRefGEM <- function(ref.data,idname){
  
  # Sub-setting the matrix by row name to extract gene short-names 
  # Pick two IDs. he gene shortnames are the same regardless of the chosen IDs. 
  # So why do these matter?
  ref1.id = row.names(subset(pData(ref.data), 
                             id == 4))
  
  ref2.id = row.names(subset(pData(ref.data), 
                             id == 5)) 
  
  # Filtering it so that only the cells with the specified IDs are shown
  ref1.filtered = ref.data[,ref1.id]
  ref2.filtered = ref.data[,ref2.id]
  
  # Create sparse GEMs for each of the references
  ref1.gem <- exprs(ref1.filtered)
  row.names(ref1.gem) <- fData(ref1.filtered)$gene_short_name
  
  ref2.gem <- exprs(ref2.filtered)
  row.names(ref2.gem) <- fData(ref2.filtered)$gene_short_name 


  # Combining the two reference GEMs into one GEM
  comb.ref.gem = cbind2(ref1.gem,ref2.gem)
  row.names(comb.ref.gem) <- fData(ref.data)$gene_short_name
  
  
  # Transform data into a format that can be used for references
  ref.transf <- SummarizedExperiment::SummarizedExperiment(
    list(counts=comb.ref.gem),
    colData=S4Vectors::DataFrame(label.main=c(ref1.filtered@phenoData@data$Main_trajectory,
                                              ref2.filtered@phenoData@data$Main_trajectory),
                                 label.sub=c(ref1.filtered@phenoData@data$Sub_trajectory_name,
                                             ref2.filtered@phenoData@data$Sub_trajectory_name))
  )
  ref.transf <- scater::logNormCounts(ref.transf)

  return(ref.transf)
  
}

IntegrateRefSample <- function(gem.seurat.obj, ref.seurat.obj){
  
  integration.list<- list(gem.seurat.obj, ref.seurat.obj) 
  
  for (i in 1:length(integration.list)) {
    integration.list[[i]] <- SCTransform(integration.list[[i]], verbose = FALSE)
  }
  
  integration.features <- SelectIntegrationFeatures(object.list = integration.list,
                                                    nfeatures = 3000)
  options(future.globals.maxSize= 3791289600)
  integration.list<- PrepSCTIntegration(object.list = integration.list, 
                                        anchor.features = integration.features, 
                                        verbose = FALSE)
  
  integration.anchors <- FindIntegrationAnchors(object.list = integration.list, 
                                                normalization.method = "SCT", 
                                                anchor.features = integration.features,
                                                verbose = FALSE)
  
  immune.combined <- IntegrateData(anchorset = integration.anchors, 
                                   normalization.method = "SCT", 
                                   verbose = FALSE)
  
  sample.ref.integrated <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
  
  return(sample.ref.integrated)
}

ClusterIntegratedData <- function(sample.ref.integrated){
  
  # t-SNE and Clustering
  sample.ref.integrated <- RunUMAP(sample.ref.integrated, reduction = "pca", dims = 1:20)
  sample.ref.integrated <- FindNeighbors(sample.ref.integrated, reduction = "pca", dims = 1:20)
  sample.ref.clustered <- FindClusters(sample.ref.integrated, resolution = 0.8)
  sample.ref.clustered$annotation = cell.annotation
  
  return(sample.ref.clustered)
  
}

CellsUnder10 <- function(cell.annotation.table.filename){
  
  cell.annotation.table <- read.table("Cell_annotation.txt",sep="\t")
  
  filtered.cells.10 <- cell.annotation.table %>% filter(Freq < 10)
  filtered.cells.10 <- filtered.cells.10[,1, drop=FALSE]
  print(filtered.cells.10, row.names = FALSE)
  
}



