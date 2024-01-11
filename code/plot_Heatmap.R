#####################
## plot_Heatmap.R  ##
#####################

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

## if you want to change conditions, you can do it here
## for example:
## gse$condition <- gsub("sg.*_","",gse$condition)

gse <- gse[,c(1,2,15,16)] ## select samples here (which two you are comparing)
gse$condition <- as.factor(gse$condition)
ddsTxi <- DESeqDataSet(gse, design = ~ condition)

mode <- "sig" ## change your parameter here

#' ! TO BE IMPLEMENTED
### save gene names (match order) ########################
if(mode == "sig") {
  sig.padj <- 0.05 ## change threshold here
  diff_res <- read.csv(file = "/Users/wanluchen/Documents/RNA_example_project/results/differential_analysis/diffRes_sg1SMARCB1_sgNTControl.csv") ## change file here
  hvg <- diff_res[which(diff_res$padj < sig.padj), "gene_id"]

  dds <- ddsTxi[hvg,]
  vsd <- vst(dds, blind=FALSE)
  heatmat <- assay(vsd)
  heatmat <- t(apply(heatmat, 1, function(x) (x-mean(x))/sd(x)))
  ## change names here if you want a different name
  ## for example:
  ## colnames(heatmat) <- c("SMARCB1", "SMARCB1", "CONTROL", "CONTROL")
  
  cl = ifelse(diff_res[which(diff_res$padj < sig.padj), "log2FoldChange"] > 0, 1, 2)
} else if(mode == "top") {
  topN <- 5000 ## change threshold here
  diff_res <- read.csv(file = "/Users/wanluchen/Documents/RNA_example_project/results/differential_analysis/diffRes_sg1SMARCB1_sgNTControl.csv") ## change file here
  diff_res <- diff_res[order(diff_res$padj),]
  hvg <- diff_res[1:topN, "gene_id"]

  dds <- ddsTxi[hvg,]
  vsd <- vst(dds, blind=FALSE)
  heatmat <- assay(vsd)
  heatmat <- t(apply(heatmat, 1, function(x) (x-mean(x))/sd(x)))
  ## change names here if you want a different name
  ## for example:
  ## colnames(heatmat) <- c("SMARCB1", "SMARCB1", "CONTROL", "CONTROL")

  cl = ifelse(diff_res[1:topN, "log2FoldChange"] > 0, 1, 2)
} else {stop("Wrong parameter!")}

# Heatmap ----------------------------------------------------------------------
ha.top = HeatmapAnnotation(
  df = data.frame(Condition = gse$condition),
  col = list(Condition = c("sgNT_Control" = "#787975", "sg1_SMARCB1" = "#721D20")), # change names and colors here
  simple_anno_size = unit(2, "mm"),
  annotation_name_gp= gpar(fontsize = 10)
)

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("#283974", "#54BCDC", "#C5C7C5", "#E1623B", "#BA3C36")) # change colors here
p <- Heatmap(heatmat, col = col_fun, name = "Significant\nGenes", top_annotation = ha.top, # change name here: "Significant\nGenes" or "Top 5k\nGenes"
             row_split = cl, cluster_row_slices = FALSE, show_row_names = FALSE,
             show_column_names = TRUE, show_row_dend = TRUE,
             show_column_dend = FALSE, column_order = c(3,4,1,2), ## change column order here, the default is "column_order = 1:ncol(heatmat)"
             ## set show_column_dend = TRUE and delete column_order to cluster columns
             row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
             column_title = "KMT2D RNA-seq", row_title = c(glue("Upregulated genes (n={length(which(cl == 1))})"), glue("Downregulated genes (n={length(which(cl == 2))})"))) # change names and numbers here
p

# save plot in pdf
dir.create("/Users/wanluchen/Documents/RNA_example_project/results/Heatmap/",
           showWarnings = FALSE, recursive = TRUE) ## change path here
pdf("/Users/wanluchen/Documents/RNA_example_project/results/Heatmap/sg1SMARCB1_sgNTControl_pval0.05.pdf", width = 5, height = 10) ## change path here
p
dev.off()

