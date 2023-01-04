
# Misc data analysis scripts

# Gene to EnsembleID conversion -------------------------------------------
library(mygene)

# create top 50 umi table
ListTopUMI <- function(gem.filtered){
  
  umi.sum <- colSums(gem.filtered[,2:ncol(gem.filtered)]) #disp number of cols
  top.umi.sum <- sort(umi.sum, decreasing = TRUE) #show how many genes expressed per col
  top50.umi <- top.umi.sum[1:50] #restrict to only show first 50 rows
  
  write.table(top50.umi, "top50_umi.tsv",
              sep = "\t", quote = F)
}

#top50.umi.list <- ListTopUMI(gem.filtered)

# convert gene list to ensembleIDs
GeneToEnsID <- function(genelist.filepath){
  
  # if there is no genelist file, make "mygenes" be this specific list of genes
  if (is.null(genelist.filepath)){
    
    genelist=c("Afp", "Trf", "Hbb.bs", "A1b", "Slc4a1", "Mt1")
  }
  
  else {
    genelist <- read.table(genelist.filepath,
                           sep = "\t",
                           row.names = NULL)[c(1)]
  }
  
  mygenes = c(genelist)
  ensemble_list <- queryMany(mygenes,
                             scopes = "symbol",
                             fields = "ensembl.gene",
                             species = "mouse")
  ensemble_list <- ensemble_list[ , c(1,5)]
  print(ensemble_list)
  
  write.table(ensemble_list, "gene_to_ensIDs.txt",
              sep = "\t", quote = F)
  
}

# create ensembleID tsv file
#ensemblID.list <- GeneToEnsID("top50_umi.tsv") # top 50 UMI list 
#ensemblID.list <- GeneToEnsID(genelist.filepath = NULL) #specific genes


# GO analysis DotPlots ----------------------------------------------------

library(tidyr)
library(ggplot2)
library(ggpubr)

# Read in GO summary table and create dotplot
GoDotplot <- function(go.summary.filepath){
  
  go.summary <- read.table(go.summary.filepath, header = FALSE, sep = "\t")
  go.summary$V6 <- as.factor(go.summary$V6)
  go.summary$V1 <- factor(go.summary$V1)
  
  # generate plot
  dotplot <- ggplot(go.summary, aes(x = V6, y = V1 , color = V5, size = V3)) + 
    geom_point() +xlab("Cluster") + ylab("Biological Process") + 
    labs(colour = "P.adjust", size = "GeneRatio") +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(size=12,colour ="black"),
          axis.title=element_text(size=14,face="bold",colour ="black"),
          legend.text=element_text(size=12),
          legend.title = element_text(colour="black", size=14, face="bold"),
          axis.line = element_line(colour = "black", size =1),
          panel.grid.major = element_line(colour = "lightgray", size =1.5), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  dotplot = dotplot+scale_color_gradient(low = "red2",  
                                         high = "mediumblue", 
                                         space = "Lab", 
                                         limit = c(0.000000000000000000000000000000000000000000000004, 0.03))
  dotplot+scale_size(range = c(2, 8))
  
}

