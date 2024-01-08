###############################
## plot_ClusteringHeatmap.R  ##
###############################

library(GenomicAlignments)
library(data.table)
library(dplyr)
library(scater)
library(GenomicRanges)
library(glue)
library(magrittr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(DESeq2)

gse <- readRDS("/Users/wanluchen/Documents/RNA_example_project/data/gse_object.rds")

## if you want to change conditions or add features, you can do it here
## for example:
gse$target <- gsub("sg.*_","",gse$condition)
gse$guide <- gsub("_.*","",gse$condition)


gse$condition <- as.factor(gse$condition)
ddsTxi <- DESeqDataSet(gse, design = ~ condition)

dds <- ddsTxi[rowSums(counts(ddsTxi) >= 10) > 2,]
vsd <- vst(dds, blind=FALSE)

heatmat <- assay(vsd)
heatmat <- t(apply(heatmat, 1, function(x) (x-mean(x))/sd(x)))

set.seed(12345)

n_group = 3
cl = kmeans(heatmat, centers = n_group, iter.max = 200)$cluster

# Heatmap ----------------------------------------------------------------------

## change annotations here
ha.top = HeatmapAnnotation(
  Target = gse$target,
  Guide = gse$guide,
  col = list(Target = c("Control" = "#33383B", "SMARCB1" = "#BBCED9", "SMARCE1" = "#dcbeff"),
             Guide = c("sg1" = "#800000", "sg2" = "#9A6324", "sg3" = "#808000", "sgGFP" = "#469990", "sgNT" = "#000075")),
  simple_anno_size = unit(4, "mm"),
  annotation_name_gp = gpar(fontsize = 10)
)

col_fun = colorRamp2(c(-2, -1.5, -1, 0, 1, 1.5, 2), c("#122240", "#1B335F", "#66A6CC", "white", "#E48367", "#670921", "#450616"))

ht_opt$message = FALSE
p <- Heatmap(heatmat, col = col_fun, name = "All\nGenes", top_annotation = ha.top,
             row_split = cl, cluster_row_slices = FALSE, show_row_names = FALSE,
             show_column_names = TRUE, show_row_dend = FALSE,
             show_column_dend = TRUE, # column_order = 1:ncol(heatmat),
             column_title = "xxx RNA samples",
             row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))

dir.create("/Users/wanluchen/Documents/RNA_example_project/results/Heatmap/",
           showWarnings = FALSE, recursive = TRUE) ## change path here
pdf("/Users/wanluchen/Documents/RNA_example_project/results/Heatmap/3clusters.pdf", width = 7, height = 9) ## change path here
p
dev.off()
