
# #!/path/to/Rscript
# library('getopt')
# 
# spec = matrix(c(
#   'verbose', 'v', 2, "integer",
#   'help', 'h', 0, "logical",
#   'gemlist', 'g', 1, "character",
#   'namelist', 'n', 1, "character",
#   'ynumber', 'y', 1, "character",
#   'runlist', 'r', 1, "character"
# ), byrow=TRUE, ncol=4)
# opt = getopt(spec)


# -------------------------------------------------------------------------

# objective of this script is to create violin plots for UMI and gene count across all samples

library(Seurat)
library(ggplot2)
library(patchwork)
library(rlang)
library(readxl)


# Read in the desired filtered gems
input.list = read_xlsx("violin_input.xlsx")
file.names <- c(input.list$path) 

# Names to be used in plotting
id.names <- c(input.list$RunID)


# Make list of gem dataframes 
runs <- list()
for (i in 1:length(id.names)){ 
  runs[[i]] <- read.table(sprintf(file.names[[i]]),
                          sep = "\t",
                          header = TRUE,
                          row.names = 1)
}


# Transpose gems, make seurat object, add field 'Name'
seurat.runs = list()
for (i in 1:length(runs)){
  runs[[i]] <- t(runs[[i]])
  seurat.runs[[i]] = CreateSeuratObject(runs[[i]], min.cells=10, project=id.names[i])
  seurat.runs[[i]][["Name"]] = id.names[i]
}

# Combine the seurat objects
run.list = c(seurat.runs[[2]],
             seurat.runs[[3]],
             seurat.runs[[4]],
             seurat.runs[[5]],
             seurat.runs[[6]],
             seurat.runs[[7]])
combined <- merge(x = seurat.runs[[1]],
                  y = run.list,
                  add.cell.ids = id.names,
                  project = "combined")


# ######
# run.list2 = c(input.list$run_list)
# combined <- merge(x = seurat.runs[[1]],
#                   y = run.list,
#                   add.cell.ids = id.names,
#                   project = "combined")


# run.list = c(input.list$run_list)
# run.list.x = run.list[1]
# run.list.y = tail(run.list,-1)
# combined <- merge(x = run.list.x,
#                   y = run.list.y,
#                   add.cell.ids = id.names,
#                   project = "combined")








# Violin plots------------------------------------------------------------------
genes = VlnPlot(object = combined, features = c("nFeature_RNA"),
                pt.size = 0,
                combine = FALSE,
                ncol = 1)

genes.plot = genes[[1]] +
  labs(y = "# genes per tixel") +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  geom_boxplot(fill = "white", width = 0.1, outlier.size=0)

genes.plot


# change scale
umi = VlnPlot(object = combined, features = c("nCount_RNA"),
              pt.size = 0,
              combine = FALSE,
              ncol = 1)

umi.plot = umi[[1]] + labs(y = "# UMIs per tixel") +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  geom_boxplot(fill='white', width = 0.1, outlier.size=0) +
  scale_y_continuous(limits = c(0, 20000)) 

umi.plot

pdf("violin_comparison_TEST.pdf", width=8.6, height=5)
wrap_plots(genes.plot, umi.plot)
dev.off()





