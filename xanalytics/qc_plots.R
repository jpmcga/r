library(ggplot2)
library(plyr)
library(magrittr)
library(tidyr)
library(OpenImageR)
library(ggpubr)
library(grid)
library(dplyr)


GemSplitXandY <- function(gem){
  # Make seperate columns for the x and y coordinates.
  gem.split <- gem %>% separate(X,
                                c("x_cord", "y_cord"),
                                sep = "x")
  return(gem.split)
}


PlotCountHistogram <- function(count, x.max, x.label = NULL){
  
  df <- data.frame(number = 1, c = count)
  
  plot <- ggplot(df, aes(x = c),
                 color ='blue',
                 xlab = x.label) + 
    geom_histogram(aes(y=..density..),
                   binwidth = x.max/20,
                   color = "black",
                   fill = "white",
                   size = 1) + 
    geom_density(alpha = .2,
                 fill = "#FF6666",
                 size = 1,
                 color = "red") +
    scale_x_continuous(name = x.label,
                       limits = c(0, x.max)) + 
    scale_y_continuous(name = "Density",
                       expand = c(0, 0)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(colour="black", size=20),
          axis.title=element_text(colour="black", size=25,face="bold"),
          legend.text=element_text(colour="black", size=20),
          legend.title = element_text(colour="black", size=20, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(colour = 'black',
                                     size=0.5,
                                     linetype='solid'),
          axis.line.y = element_line(colour = 'black',
                                     size=0.5,
                                     linetype='solid'))
  return(plot)
}


PlotHeatMaps <- function(gem,
                         count,
                         max = 20000,
                         tixel.size = 3,
                         name = NULL,
                         image.file = NULL){
  
  gem <- GemSplitXandY(gem)
  
  if (! is.null(image.file)){
    image <- OpenImageR::readImage(image.file)
    image.raster <- rasterGrob(image, width = unit(1,"npc"),
                               height=unit(1,"npc"),
                               interpolate = FALSE)
  }
  
  plot <- ggplot(gem,
                 aes(x = as.numeric(x_cord), 
                     y = as.numeric(y_cord),
                     color = count)) +
    scale_color_gradientn(colours = c("blue","green", "red"),
                          limits = c(0, max),
                          oob = scales::squish) +
    ggtitle(sprintf("Heatmap of %s Counts", name)) +
    {if (! is.null(image.file))
      annotation_custom(image.raster,
                        xmin=-Inf,
                        xmax=Inf,
                        ymin=-Inf,
                        ymax=Inf)
    } +
    guides(colour = guide_colourbar(barwidth = 1,
                                    barheight = 30)) +
    geom_point(shape = 15, size = tixel.size) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X",
                       limits = c(NA, NA),
                       expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y",
                    limits = c(NA, NA),
                    expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51), ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face = "bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  return(plot)
  
}



