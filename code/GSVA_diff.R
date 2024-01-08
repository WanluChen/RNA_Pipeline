##################
## GSVA_diff.R  ##
##################

# Usage:
# gse_object: gse object obtained in process_data.R
# gene_set, sub_set: pathways; if no subset then don't include it
# Note that there is a helper function to show the available collections:
# msigdbr::msigdbr_collections() %>% as.data.frame()
# file_name: output file name (you will get two files: .csv and .pdf)
# show_name: whether showing the pathway(gene set) names
# rowname_size: font size of pathway(gene set) names, only used when show_name = T
# diff_analysis: whether perform differential analysis
# adj_pval: adjusted p values for differentially expressed pathways, only used when diff_analysis = T
# sample_of_interest: the index of samples of interest; perform differential analysis on [samples of interest] vs [all remaining samples], only used when diff_analysis = T

library(GSVA)
library(magrittr)
library(org.Hs.eg.db) ## human
# if (!require("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db");library(org.Mm.eg.db) ## mouse
library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(limma)

GSVA_diff <- function(gse_object, gene_set, sub_set = NULL, file_name, show_name, rowname_size = 5, diff_analysis = T, adj_pval = 0.05, sample_of_interest){
  
  ddsTxi <- DESeqDataSet(gse_object, design = ~ 1)
  dds <- ddsTxi[rowSums(counts(ddsTxi) >= 10) > 2,]
  dds_norm <- vst(dds, blind=FALSE)
  vst_df <- assay(dds_norm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ensembl_id")
  
  gene_sets <- msigdbr::msigdbr(
    species = "Homo sapiens",
    category = gene_set,
    subcategory = sub_set)
  pathway_list <- split(
    gene_sets$entrez_gene,
    gene_sets$gs_name
  )
  
  mapped_df <- data.frame(
    "entrez_id" = mapIds(
      org.Hs.eg.db,
      keys = vst_df$ensembl_id,
      keytype = "ENSEMBL",
      column = "ENTREZID",
      multiVals = "first"
    )
  ) %>%
    dplyr::filter(!is.na(entrez_id)) %>%
    tibble::rownames_to_column("Ensembl") %>%
    dplyr::inner_join(vst_df, by = c("Ensembl" = "ensembl_id"))
  gene_means <- rowMeans(mapped_df %>% dplyr::select(-Ensembl, -entrez_id))
  mapped_df <- mapped_df %>%
    dplyr::mutate(gene_means) %>%
    dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())
  filtered_mapped_df <- mapped_df %>%
    dplyr::arrange(dplyr::desc(gene_means)) %>%
    dplyr::distinct(entrez_id, .keep_all = TRUE)
  filtered_mapped_matrix <- filtered_mapped_df %>%
    dplyr::select(-Ensembl, -gene_means) %>%
    tibble::column_to_rownames("entrez_id") %>%
    as.matrix()
  
  gsvaPar <- gsvaParam(
    exprData = filtered_mapped_matrix,
    geneSets = pathway_list,
    minSize = 10,
    maxSize = 500,
    kcdf = "Gaussian",
    maxDiff = TRUE,
    absRanking = FALSE)
  gsva_results <- gsva(gsvaPar, verbose=FALSE)
  
  tmp <- rep("Ref", ncol(gse_object))
  tmp[sample_of_interest] <- "Focus"
  gse_object$condition <- relevel(as.factor(tmp), ref = "Ref")

  mod <- model.matrix(~ factor(gse_object$condition))
  fit <- lmFit(gsva_results, mod)
  fit <- eBayes(fit)
  res <- decideTests(fit, p.value = adj_pval)
  
  tt <- topTable(fit, coef=2, n=Inf) # control mutant
  DEpwys <- rownames(tt)[tt$adj.P.Val <= adj_pval]
  pdf(paste0(output_path, file_name, "_diff_volcano.pdf"), width = 5, height = 5)
  plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
       main="", xlab="GSVA enrichment score difference",
       ylab=expression(-log[10]~~Raw~P-value))
  abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= adj_pval])),
         col=grey(0.5), lwd=1, lty=2)
  points(tt$logFC[match(DEpwys, rownames(tt))],
         -log10(tt$P.Value[match(DEpwys, rownames(tt))]),
         pch=".", cex=5, col="darkred")
  text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= adj_pval])),
       paste0(adj_pval*100, "% FDR"), pos=3)
  dev.off()
  
  pathway_list <- pathway_list[DEpwys]
  gsvaPar <- gsvaParam(
    exprData = filtered_mapped_matrix,
    geneSets = pathway_list,
    minSize = 10,
    maxSize = 500,
    kcdf = "Gaussian",
    maxDiff = TRUE,
    absRanking = FALSE)
  gsva_results <- gsva(gsvaPar, verbose=FALSE)
  
  ht_opt$message = FALSE
  colvec <- setNames(col_pool[1:length(unique(gse_object$condition))], unique(gse_object$condition))
  ha.top = HeatmapAnnotation(
    Condition = gse_object$condition,
    col = list(Condition = colvec),
    simple_anno_size = unit(4, "mm"),
    annotation_name_gp = gpar(fontsize = 10)
  )
  col_fun = colorRamp2(c(min(gsva_results), (min(gsva_results) + max(gsva_results))/2, max(gsva_results)), c("#65A67F", "white", "#D42D00"))
  lgd = Legend(col_fun = col_fun, title = "GSVA\nScore")
  ht_list = Heatmap(gsva_results, col = col_fun, name = "GSVA\nScore", top_annotation = ha.top,
                    show_row_names = show_name, show_column_names = TRUE,
                    show_row_dend = TRUE, column_order = c(sample_of_interest, setdiff(1:ncol(gsva_results), sample_of_interest)),
                    column_title = file_name,
                    row_names_gp = gpar(fontsize = rowname_size), column_names_gp = gpar(fontsize = 10))
  
  pdf(paste0(output_path, file_name, "_diff.pdf"), width = 7, height = 10) ## change plot size here
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "left",
       annotation_legend_side = "left", padding = unit(c(2, 2, 2, 15), "mm"))
  dev.off()
  
  return(res)
  
}

col_pool <- c("#A50226", "#D73028", "#F56D43", "#FEAF61", "#FEE091", "#FFFFBE", "#DFF3F9", "#ACD8E9", "#75ACD2", "#4774B5", "#303794")

output_path <- "/Users/wanluchen/Documents/RNA_example_project/results/GSVA/" ## change path here
gse <- readRDS("/Users/wanluchen/Documents/RNA_example_project/data/gse_object.rds") ## change path here
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)


res <- GSVA_diff(gse_object = gse, gene_set = "H", file_name = "Hallmark.Control_REST",
                 show_name = T, rowname_size = 10, diff_analysis = T, adj_pval = 0.05, sample_of_interest = grep("Control", gse$condition))
summary(res) # check numbers of up/down-regulated pathways
