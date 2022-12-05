library(dplyr)

dir <- "/Users/jamesmcgann/Atlasxomics/Data/AtlasXAnalytics/experiments/D04-jm"  
setwd(dir)

gem <- read.table(file = 'text_files/D04_stdata_names.tsv',
                  sep = '\t',
                  header = TRUE,
                  stringsAsFactors=FALSE)

col1 <- gem$X
unique.cols <- which(colSums(gem) == 1)
minus1_fil <- gem[, colSums(gem[-1]) > 10]
final = data.frame(X = gem$X, gem[, colSums(gem[-1]) > 10])

new <- gem[c(1,-1)]

gem3 <- merge(gem$X, gem[-1])