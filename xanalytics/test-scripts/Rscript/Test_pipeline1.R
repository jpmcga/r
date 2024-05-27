
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
  'gem', 'g', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)


# -------------------------------------------------------------------------


library(Seurat)

source("../DBiT datasets Seurat/preprocessing.R")
source("../DBiT datasets Seurat/qc_plots.R")


IMAGE.FILE <- opt$image


# Preprocessing ---------------------------------------------------------------

# Filter Gene Expression Matrix (GEM) to include only tixels on tissue
gem.filtered <- FilterTixelsOnTissue(opt$position,
                                     opt$gem)

write.table(gem.filtered,
            "test_filtered_gem.tsv",
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)


# -------------------------------------------------------------------------


# Rscript --vanilla Test_pipeline1.R --image=test_ROI.png --BWimage=test_ROI_BW.png --position=test_position.txt --gem=test_stdata_names.tsv 






