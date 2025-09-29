library(ComplexHeatmap)
library(circlize)
library(magick)
library(DESeq2)
library(ggplot2)
library(glue)

setwd("/Volumes/KYOYA/Genomics/ToskaLab/CaiLab/FusionTFE3_RNAseq")
dir.create("results/Heatmap/", showWarnings = FALSE, recursive = TRUE)
gse <- readRDS("data/gse_object.rds")

plotHeatmap <- function(gse, n_group, file_name, cl.order = NULL){
  gse$condition <- as.factor(gse$condition)
  ddsTxi <- DESeqDataSet(gse, design = ~ condition)
  dds <- ddsTxi[rowSums(counts(ddsTxi) >= 10) > 2,]
  vsd <- vst(dds, blind=FALSE)
  heatmat <- assay(vsd)
  heatmat <- t(apply(heatmat, 1, function(x) (x-mean(x))/sd(x)))
  
  set.seed(12345)
  cl = kmeans(heatmat, centers = n_group, iter.max = 200)$cluster
  
  # reorder clusters
  if(!is.null(cl.order)) {
    k = 1
    for (i in cl.order) {
      cl[which(cl==i)] <- paste0("cl_", k)
      k <- k+1
    }
  } else {
    tmp <- names(cl)
    cl <- paste0("cl_", cl)
    names(cl) <- tmp
  }
  
  # # save genes by cluster
  # for (i in 1:n_group) {
  #   cl.df <- rowData(gse)[names(cl[which(cl == paste0("cl_", i))]),] %>% as.data.frame()
  #   cl.df$entrezid <- as.character(cl.df$entrezid)
  #   cl.df$tx_ids <- as.character(cl.df$tx_ids)
  #   write.csv(cl.df, glue("results/Heatmap/{file_name}_{n_group}clusters_cl{i}.csv"), row.names = F)
  # }
  
  # add number of genes
  tmp <- names(cl)
  cl <- glue("{cl}(n={table(cl)[cl]})")
  names(cl) <- tmp

  # Heatmap -------------
  ha.top = HeatmapAnnotation(
    Condition = gse$condition,
    col = list(Condition = c("KO" = "#000000", "NONO_TFE3" = "#B10026", "NONOdCCD_TFE3" = "#fd8d3c", "SFPQ_TFE3" = "#0d4a70",
                             "SFPQdCCD_TFE3" = "#6cb0d6", "NONOccmut_TFE3" = "#feb24c", "CM" = "#a0a0a4")),
    simple_anno_size = unit(4, "mm"),
    annotation_name_gp = gpar(fontsize = 10)
  )
  col_fun = colorRamp2(c(-2, -1.5, -1, 0, 1, 1.5, 2), c("#122240", "#1B335F", "#66A6CC", "white", "#E48367", "#670921", "#450616"))
  ht_opt$message = FALSE
  p <- Heatmap(heatmat, col = col_fun, name = "All\nGenes", top_annotation = ha.top,
               row_split = cl, cluster_row_slices = FALSE, show_row_names = FALSE,
               show_column_names = TRUE, show_row_dend = FALSE,
               show_column_dend = TRUE, # column_order = 1:ncol(heatmat),
               column_title = file_name,
               row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
  pdf(glue("results/Heatmap/{file_name}_{n_group}clusters.pdf"), width = 7, height = 8)
  print(p)
  dev.off()
}

# no EBSS
gse1 <- gse[,-which(gse$condition == "EBSS")]
plotHeatmap(gse1, n_group = 5, file_name = "noEBSS")
plotHeatmap(gse1, n_group = 6, file_name = "noEBSS")
plotHeatmap(gse1, n_group = 7, file_name = "noEBSS")

## get genes by cluster
files <- list.files("results/DEseq2", pattern = ".csv", full.names = T)
files <- files[-grep("EBSS", files, value = F)]
refgene <- read.csv("results/Heatmap/noEBSS_6clusters_cl3.csv")
res <- refgene
for (f in files) {
  tmp <- read.csv(f, row.names = "gene_id")
  tmp$X <- NULL
  tmp <- tmp[refgene$gene_id, c("log2FoldChange", "padj")]
  colnames(tmp) <- paste0(gsub(".csv", "", gsub("results/DEseq2/", "", f)), ".", c("log2FC", "padj"))
  res <- cbind(res, tmp)
}
gse1 <- gse[refgene$gene_id, -which(gse$condition == "EBSS")]
res <- cbind(res,  assays(gse1)[[2]])
write.csv(res, "results/Heatmap/noEBSS_6clusters_cl3.csv")


# no EBSS/KO
gse1 <- gse[,-which(gse$condition == "EBSS" | gse$condition == "KO")]
plotHeatmap(gse1, n_group = 5, file_name = "noEBSS_noKO")
plotHeatmap(gse1, n_group = 6, file_name = "noEBSS_noKO")
plotHeatmap(gse1, n_group = 7, file_name = "noEBSS_noKO")


# (k = 3/4/5) with KO, 118CM, 106 NONOTFE3 and 112 SFPQTFE3
gse1 <- gse[,which(gse$condition == "KO" | gse$condition == "CM" | gse$condition == "NONO_TFE3" | gse$condition == "SFPQ_TFE3")]
plotHeatmap(gse1, n_group = 3, file_name = "KO_CM_NONO_SFPQ")
plotHeatmap(gse1, n_group = 4, file_name = "KO_CM_NONO_SFPQ")
plotHeatmap(gse1, n_group = 5, file_name = "KO_CM_NONO_SFPQ")


# (k = 4/5/6) with 118CM, 106 NONOTFE3, 112 SFPQTFE3, 111 NONOdCCD and 116 SFPQdCCD
gse1 <- gse[,which(gse$condition == "CM" | gse$condition == "NONO_TFE3" | gse$condition == "SFPQ_TFE3" |
                    gse$condition == "NONOdCCD_TFE3" | gse$condition == "SFPQdCCD_TFE3")]
plotHeatmap(gse1, n_group = 4, file_name = "CM_NONO_SFPQ_dCCD")
plotHeatmap(gse1, n_group = 5, file_name = "CM_NONO_SFPQ_dCCD")
plotHeatmap(gse1, n_group = 6, file_name = "CM_NONO_SFPQ_dCCD")


# (k = 3/4/5) with 106 NONOTFE3, 111NONOdCCD and 119 NONOccmut
gse1 <- gse[,which(gse$condition == "NONO_TFE3" | gse$condition == "NONOdCCD_TFE3" | gse$condition == "NONOccmut_TFE3")]
plotHeatmap(gse1, n_group = 3, file_name = "NONO_dCCD_ccmut")
plotHeatmap(gse1, n_group = 4, file_name = "NONO_dCCD_ccmut")
plotHeatmap(gse1, n_group = 5, file_name = "NONO_dCCD_ccmut")

# version 2 --------------------------------------------------------------------
setwd("/Volumes/KYOYA/Genomics/ToskaLab/CaiLab/FusionTFE3_RNAseq")
dir.create("results/Heatmap/", showWarnings = FALSE, recursive = TRUE)

plotHeatmap <- function(gse, n_group, file_name, cl.order = NULL){
  gse$condition <- as.factor(gse$condition)
  ddsTxi <- DESeqDataSet(gse, design = ~ condition)
  dds <- ddsTxi[rowSums(counts(ddsTxi) >= 10) > 2,]
  if(nrow(dds) < 1000) {
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  } else {
    vsd <- vst(dds, blind=FALSE)
  }
  heatmat <- assay(vsd)
  heatmat <- t(apply(heatmat, 1, function(x) (x-mean(x))/sd(x)))
  
  set.seed(12345)
  cl = kmeans(heatmat, centers = n_group, iter.max = 200)$cluster
  # reorder clusters
  if(!is.null(cl.order)) {
    k = 1
    for (i in cl.order) {
      cl[which(cl==i)] <- paste0("cl_", k)
      k <- k+1
    }
  } else {
    tmp <- names(cl)
    cl <- paste0("cl_", cl)
    names(cl) <- tmp
  }

  # save genes by cluster
  # for (i in 1:n_group) {
  #   cl.df <- rowData(gse)[names(cl[which(cl == paste0("cl_", i))]),] %>% as.data.frame()
  #   cl.df$entrezid <- as.character(cl.df$entrezid)
  #   cl.df$tx_ids <- as.character(cl.df$tx_ids)
  #   write.csv(cl.df, glue("results/Heatmap/RNA_{file_name}_{n_group}clusters_cl{i}.csv"), row.names = F)
  # }
  
  # Heatmap -------------
  ha.top = HeatmapAnnotation(
    Condition = gse$condition,
    col = list(Condition = c("KO" = "#000000", "NONO_TFE3" = "#B10026", "NONOdCCD_TFE3" = "#fd8d3c", "SFPQ_TFE3" = "#0d4a70",
                             "SFPQdCCD_TFE3" = "#6cb0d6", "NONOccmut_TFE3" = "#feb24c", "CM" = "#a0a0a4")),
    simple_anno_size = unit(4, "mm"),
    annotation_name_gp = gpar(fontsize = 10)
  )
  col_fun = colorRamp2(c(-2, -1.5, -1, 0, 1, 1.5, 2), c("#122240", "#1B335F", "#66A6CC", "white", "#E48367", "#670921", "#450616"))
  ht_opt$message = FALSE
  tmp <- names(cl)
  cl <- glue("{cl}(n={table(cl)[cl]})")
  names(cl) <- tmp
  p <- Heatmap(heatmat, col = col_fun, name = "z-score", top_annotation = ha.top,
               row_split = cl, cluster_row_slices = FALSE, show_row_names = FALSE,
               show_column_names = TRUE, show_row_dend = FALSE,
               show_column_dend = FALSE, column_order = 1:ncol(heatmat),
               column_title = file_name,
               row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))
  pdf(glue("results/Heatmap/RNA_{file_name}_{n_group}clusters.pdf"), width = 7, height = 8)
  print(p)
  dev.off()
}

# KO, 118CM, 106 NONOTFE3 and 112 SFPQTFE3
gse1 <- readRDS("results/DEseq2/KO-else_gse.rds")
gse1 <- gse1[,c(which(gse1$condition == "KO"),
                which(gse1$condition == "CM"),
                which(gse1$condition == "NONO_TFE3"),
                which(gse1$condition == "SFPQ_TFE3"))]
plotHeatmap(gse1, n_group = 5, file_name = "KO_CM_NONO_SFPQ", cl.order = c(5,4,1,3,2))

# NONO/SFPQ dCCD, NONO/SFPQ
gse1 <- readRDS("results/DEseq2/NONO_SFPQ_dCCD-NONO_SFPQ_gse.rds")
gse1 <- gse1[,c(which(gse1$condition == "NONO_TFE3"),
                which(gse1$condition == "SFPQ_TFE3"),
                which(gse1$condition == "NONOdCCD_TFE3"),
                which(gse1$condition == "SFPQdCCD_TFE3"))]
plotHeatmap(gse1, n_group = 4, file_name = "NONO_SFPQ_dCCD", cl.order = c(4,2,3,1))
 
# NONO/SFPQ dCCD /NONOccmut, NONO/SFPQ
gse1 <- readRDS("results/DEseq2/NONO_SFPQ_dCCD_ccmut-NONO_SFPQ_gse.rds")
gse1 <- gse1[,c(which(gse1$condition == "NONO_TFE3"),
                which(gse1$condition == "SFPQ_TFE3"),
                which(gse1$condition == "NONOdCCD_TFE3"),
                which(gse1$condition == "SFPQdCCD_TFE3"),
                which(gse1$condition == "NONOccmut_TFE3"))]
plotHeatmap(gse1, n_group = 4, file_name = "NONO_SFPQ_dCCD_ccmut", cl.order = c(2,4,1,3))

# NONO dCCD, ccmut
gse1 <- readRDS("results/DEseq2/NONO-dCCD_ccmut_gse.rds")
gse1 <- gse1[,c(which(gse1$condition == "NONO_TFE3"),
                which(gse1$condition == "NONOdCCD_TFE3"),
                which(gse1$condition == "NONOccmut_TFE3"))]
plotHeatmap(gse1, n_group = 2, file_name = "NONO_dCCD_ccmut")
 


