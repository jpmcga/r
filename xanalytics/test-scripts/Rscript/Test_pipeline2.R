
#!/path/to/Rscript
library('getopt')
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).

spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help', 'h', 0, "logical",
  'image', 'i', 1, "character",
  'BWimage', 'b', 1, "character",
  'filteredgem', 'f', 1, "character",
  'mtpattern', 'm', 1, "character",
  'genelist', 'l', 1, "character",
  'id', 'd', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# -------------------------------------------------------------------------

library(Seurat)
source("../DBiT datasets Seurat/preprocessing.R")
source("../DBiT datasets Seurat/qc_plots.R")
source("../DBiT datasets Seurat/processing.R")
source("../DBiT datasets Seurat/clustering_plots.R")
source("../DBiT datasets Seurat/gene_plots.R")

run.id <- opt$id
IMAGE.FILE <- opt$image

# Remove genes with < 10 tixels/per from GEM
gem.filtered.gt10 <- FilterLtNCells(opt$filteredgem, min.tixels = 10)

# Get number of UMI counts in each tixel
umi.count <- GetUMIcount(gem.filtered.gt10)

# Get number of gene counts in each tixel
gene.count <- GetGeneCount(gem.filtered.gt10)

# Get number of tixels on tissue
num.tixels <- nrow(gem.filtered.gt10)

# Print stats
print("Number of Tixels on Tissue")
print(num.tixels)
avg.umi <- mean(umi.count)
avg.genes <- mean(gene.count) 
print("Average UMI per tixel")
print(avg.umi)
print("Average genes per tixel")
print(avg.genes)

# Create stats txt file
file.name <- paste(run.id, "run_info.txt", sep = "_")
run.stats <- file(file.name)
writeLines(c("Number of tixels on tissue:",
             num.tixels,
             "Average UMI per tixel:", 
             avg.umi,
             "Average genes per tixel:",
             avg.genes), run.stats)
close(run.stats)


#  Processing ------------------------------------------------------------------

# Create Seurat object and tixel cluster df
print("Making Seurat Oject")
gem.seurat.obj <- GemToSeuratObj(gem.filepath = opt$filteredgem,
                                 mt.pattern = opt$mtpattern) 
# "^mt." for mouse, "^MT-" for human

print("Making tixel cluster df")
tixel.cluster.df <- GetTixelClusterDf(gem.seurat.obj)
max.clusters <- max(c(tixel.cluster.df$count)-1) # total number of clusters


# Get Differentially Expressed Genes (DEGs)
print("Getting DEGs")
degs <- GetDEGs(gem.seurat.obj)

# Get top n degs for each cluster
print("Finding top n Markers")
deg.top20 <- GetTopnDEGMarkers(degs)



# Plots -------------------------------------------------------------------
pdf.name <- paste(run.id, "report.pdf", sep = "_")

pdf(file=pdf.name, width = 8.6, height = 8.6)

# QC histograms and heatmaps 
print(PlotCountHistogram(gene.count, 8000, "Gene"))
print(PlotCountHistogram(umi.count, 30000, "UMI"))

print(PlotHeatMaps(gem.filtered.gt10,
                   count = gene.count,
                   name = "Gene",
                   max = 6000,
                   tixel.size = 1.5,
                   image.file = IMAGE.FILE))
print(PlotHeatMaps(gem.filtered.gt10,
                   count = umi.count,
                   name = "UMI",
                   max = 15000,
                   tixel.size = 1.5,
                   image.file = IMAGE.FILE))
print(PlotHeatMaps(gem.filtered.gt10,
                   count = gene.count,
                   name = "Gene",
                   max = 6000,
                   tixel.size = 1.5,
                   image.file = opt$BWimage))
print(PlotHeatMaps(gem.filtered.gt10,
                   count = umi.count,
                   name = "UMI",
                   max = 15000,
                   tixel.size = 1.5,
                   image.file = opt$BWimage))

# Clustering plots 
print(PlotUMAPClusters(gem.seurat.obj,
                       title = "UMAP"))
print(PlotSpatialClusters(tixel.cluster.df,
                          title = "UMAP",
                          tixel.size = 2.5))
print(PlotSpatialClusters(tixel.cluster.df,
                          title = "UMAP",
                          tixel.size = 1.5,
                          image.file = opt$image))
print(PlotSpatialClusters(tixel.cluster.df,
                          title = "UMAP",
                          tixel.size = 1.5,
                          image.file = opt$BWimage))

# Individual cluster spatial heatmaps
for (i in 0:max.clusters){
  
  cluster <- tixel.cluster.df[!(tixel.cluster.df$count!=i),]
  
  print(PlotSpatialClusters(cluster,
                            title = "UMAP",
                            tixel.size = 1.5,
                            image.file = opt$BWimage))
}

# DEG plots 

print("plotting DEG figures")
PlotDEGTopnHeatMap(gem.seurat.obj, deg.top20)

#Top n DEG table
deg.table.name <- paste(run.id, "Top20DEG.xlsx", sep = "_")
write.xlsx(deg.top20, deg.table.name, sep="\t")


# write.xlsx(deg.top20, "Top20DEG.xlsx",sep="\t")


# Individual gene plots 
gem.plotting <- as.data.frame(t(gem.seurat.obj[["SCT"]]@data))

gene.list <-read.table(opt$genelist, sep = "", header = TRUE)
genes <- c(gene.list$x)

gem.plotting$tixels = row.names(gem.plotting)
gem.plotting <- gem.plotting %>% separate(tixels,
                                          c("x_cord", "y_cord"),
                                          sep = "x")

for (gene in genes){
  if (gene %in% colnames(gem.plotting)){
    print(sprintf("plotting %s", gene))
    try(print(PlotGeneUMI(gem.plotting,
                          gene,
                          image.file = opt$BWimage)))
  }
}
dev.off()


# # -------------------------------------------------------------------------


# Rscript --vanilla Test_pipeline2.R --image=test_ROI.png --BWimage=test_ROI_BW.png  --filteredgem=test_filtered_gem.tsv --mtpattern=^mt. --genelist=ind_gene_list.txt --id=Dxx






