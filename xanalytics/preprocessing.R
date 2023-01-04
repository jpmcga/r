library(plyr)
library(magrittr)
library(tidyr)
library(dplyr)


FilterTixelsOnTissue <- function(position.filepath, gem.filepath){ 
  
  # Filter out tixels in Gene Expression Matrix (GEM) that are not on tissue.
  # Takes position file from image processing script and GEM tsv,
  # outputs filtered GEM as dataframe object.
  
  # Read file listing tixels on tissue.
  tixels.on.tissue <- read.table(position.filepath,
                                 sep =",",
                                 header = FALSE,
                                 dec =".",
                                 stringsAsFactors = FALSE)
  
  # Convert tixels from dataframe to character vector, remove first value, "NA".
  tixels.on.tissue <- as.character(tixels.on.tissue[1,])[-1]
  
  # Read in gem .tsv file as dataframe
  gem <- read.table(gem.filepath,
                    sep = '\t',
                    header = TRUE,
                    stringsAsFactors = FALSE)
  
  # Filter GEM, keeping only rows with tixel coordinates found in on tissue.
  gem.filtered <- gem[gem$X %in% tixels.on.tissue,]
  
  return(gem.filtered)
}


FilterLtNCells <- function(gem, min.tixels = 10){
  # Removes genes (columns) from gem dataframe object which occur in fewer than
  # min.tixels.
  
  ##
  
  # Read in gem .tsv file as dataframe
  gem <- read.table(gem,
                    sep = '\t',
                    header = TRUE,
                    stringsAsFactors = FALSE)
  
  ##
  
  gem.filtered <- data.frame(X = gem$X,
                             gem[, colSums(gem[-1]) > min.tixels])
  
  return(gem)
} 


GetUMIcount <- function(gem){
  
  # From gem dataframe return total UMI counts for each tixel
  
  umi.count <- rowSums(gem[-1]) # First column = tixels
  
  return(umi.count)
}


GetGeneCount <- function(gem){
  
  # Remove first column (tixel positions), turn into binary (if UMI 1, else 0).
  gem.binary <- gem[-1] %>% mutate_all(as.logical)
  
  # Get total gene count for each tixel.
  gene.count <- rowSums(gem.binary)
  
  return(gene.count)
}