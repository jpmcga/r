  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(magrittr)
  library(tidyr)
  library(Matrix)
  library(grid)
  library(OpenImageR)
  
  
  MakeGEMForPlotting <- function(seuratobj = NULL,
                                 filepath = NULL,
                                 normalize = TRUE){
  
    # Input gene expression matrix (GEM) from file, output normalized GEM as
    # dataframe with rows for x and y coordinates for each tixel.
  
    # Read filtered GEM from file.
    gem <- read.table(filepath,
                      sep = '\t',
                      header = TRUE,
                      row.names = 1)
    
    if (normalize){
      # Transpose GEM (rows = genes, columns = tixels) so
      # proper input for Seurat.
      gem.transposed <- t(gem)
      experiment <- CreateSeuratObject(gem.transposed) # No filtering w min.cells
      experiment <- NormalizeData(experiment,
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000)
      # Access GEM, re-transpose (columns = genes, rows = tixels),
      # create dataframe.
      gem <- as.data.frame(t(experiment[["RNA"]]@data))
    }
  
    # Duplicate row.names column (i.e tixels) and split into two columns
    # giving x and y coordinates. 
    gem$tixels = row.names(gem)
    gem <- gem %>% separate(tixels,
                            c("x_cord", "y_cord"),
                            sep = "x")
    return(gem)
  }
  
  
  PlotGeneUMI <- function(normalized_gem, gene, image.file = NULL) {

    # Outputs heatmaps for individuals genes given a gene data table (made
    # via MakeGeneTable) and a gene name. 
      
      if (! is.null(image.file)){
        image <- OpenImageR::readImage(image.file)
        image.raster <- rasterGrob(image, width = unit(1,"npc"),
                                 height=unit(1,"npc"),
                                 interpolate = FALSE)
      }
   
    
      plot <- ggplot(normalized_gem,
              aes(x = as.numeric(x_cord), y = as.numeric(y_cord),
              color = unlist(normalized_gem[gene]))) +
        scale_color_gradientn(colours = c("black", "green", "yellow"),
                              oob = scales::squish) +
        ggtitle(sprintf("%s", gene)) +
        {if (! is.null(image.file))
          annotation_custom(image.raster,
                            xmin=-Inf,
                            xmax=Inf,
                            ymin=-Inf,
                            ymax=Inf)
        } +
        guides(colour = guide_colourbar(barwidth = 1,
                                        barheight = 30)) +
        geom_point(shape = 15, size = 2) +
        expand_limits(x = 0, y = 0) +
        scale_x_continuous(name = "X",
                          limits = c(NA, NA),
                          expand = expansion(mult = c(-0.013, -0.013))) +
        scale_y_reverse(name = "Y", limits = c(NA, NA),
                        expand = expansion(mult = c(-0.013, 0.008))) +
        coord_equal(xlim = c(0,51), ylim = c(51,1)) +
        theme(plot.title = element_text(hjust = 0.5,
                                        size = 40,
                                        face = "bold"),
              axis.text=element_text(colour = "black",
                                      size = 30),
              axis.title=element_text(colour = "black",
                                      size = 30,
                                      face = "bold"),
              legend.text=element_text(colour = "black", size = 30),
              legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(size = 1, fill = NA)) 

      return(plot)
  }
  
  
  
  