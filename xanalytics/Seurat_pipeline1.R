
#!/path/to/Rscript
library('getopt')
# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).

spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help', 'h', 0, "logical",
  'image', 'i', 1, "character",
  'BWimage', 'b', 1, "character",
  'position', 'p', 1, "character",
  'gem', 'g', 1, "character",
  'id', 'd', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# -------------------------------------------------------------------------


library(Seurat)

source("../preprocessing.R")
source("../qc_plots.R")
source("../processing.R")
source("../clustering_plots.R")
source("../gene_plots.R")

run.id <- opt$id
IMAGE.FILE <- opt$image

# Preprocessing ---------------------------------------------------------------


# Read in gem .tsv file as dataframe
gem <- read.table(opt$gem,
                  sep = '\t',
                  header = TRUE,
                  stringsAsFactors = FALSE)

# Get number of tixels and UMI/gene counts in each tixel BEFORE filtering
umi.count <- GetUMIcount(gem)
gene.count <- GetGeneCount(gem)
num.tixels <- nrow(gem)
avg.umi <- mean(umi.count)
avg.genes <- mean(gene.count) 

# Create UNFILTERED stats txt file
file.name <- paste(run.id, "unfiltered_gem_info.txt", sep = "_")
run.stats <- file(file.name)
writeLines(c("Run ID:", 
             opt$id,
             "Number of tixels on tissue before filtering:",
             num.tixels,
             "Average UMI per tixel before filtering:", 
             avg.umi,
             "Average genes per tixe before filteringl:",
             avg.genes), run.stats)
close(run.stats)


# Filter Gene Expression Matrix (GEM) to include only tixels on tissue
gem.filtered <- FilterTixelsOnTissue(opt$position,
                                     opt$gem)



filtered.matrix.name <- paste(run.id, "filtered_gem.tsv", sep = "_")
write.table(gem.filtered,
            filtered.matrix.name,
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)



# # -------------------------------------------------------------------------


# Rscript --vanilla Seurat_pipeline1.R --image=Dxx_ROI.png --BWimage=Dxx_ROI_BW.png --position=Dxx_position_list.txt --gem=Dxx_stdata_names.tsv --id=Dxx







