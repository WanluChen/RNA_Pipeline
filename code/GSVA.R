#############
## GSVA.R  ##
#############

# Usage:
# gse_object: gse object obtained in process_data.R
# gene_set, sub_set: pathways; if no subset then don't include it
  # Note that there is a helper function to show the available collections:
  # msigdbr::msigdbr_collections() %>% as.data.frame()
# file_name: output file name (you will get two files: .csv and .pdf)
# show_name: whether showing the pathway(gene set) names
# rowname_size: font size of pathway(gene set) names, only used when show_name = T
# pick_pathway: whether using user-specific pathways
# user_pathway: user-specific pathways, only used when pick_pathway = F


library(GSVA)
library(magrittr)
library(org.Hs.eg.db) ## human
# if (!require("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db");library(org.Mm.eg.db) ## mouse
library(ComplexHeatmap)
library(circlize)
library(DESeq2)

GSVA <- function(gse_object, gene_set, sub_set = NULL, file_name, show_name, rowname_size = 5, pick_pathway = F, user_pathway = NULL){
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
  if(pick_pathway) {
    pathway_list <- pathway_list[user_pathway]
  }
  
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
  
  gsva_results %>%
    as.data.frame() %>%
    tibble::rownames_to_column("pathway") %>%
    write.csv(paste0(output_path, file_name, ".csv"), row.names = F)
  rownames(gsva_results) <- gsub("HALLMARK_|BIOCARTA_|KEGG_|PID_", "", rownames(gsva_results))
  
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
                    show_row_dend = TRUE,
                    show_column_dend = TRUE, # column_order = 1:ncol(gsva_results), ## set column_order here if show_column_dend = FALSE
                    column_title = file_name,
                    row_names_gp = gpar(fontsize = rowname_size), column_names_gp = gpar(fontsize = 10))
  
  pdf(paste0(output_path, file_name, ".pdf"), width = 7, height = 10) ## change plot size here
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "left",
       annotation_legend_side = "left", padding = unit(c(2, 2, 2, 15), "mm"))
  dev.off()
  
}

col_pool <- c("#A50226", "#D73028", "#F56D43", "#FEAF61", "#FEE091", "#FFFFBE", "#DFF3F9", "#ACD8E9", "#75ACD2", "#4774B5", "#303794")

output_path <- "/Users/wanluchen/Documents/RNA_example_project/results/GSVA/" ## change path here
gse <- readRDS("/Users/wanluchen/Documents/RNA_example_project/data/gse_object.rds") ## change path here
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

## examples
GSVA(gse_object = gse, gene_set = "H", file_name = "Hallmark", show_name = T, rowname_size = 10)
GSVA(gse_object = gse, gene_set = "C2", file_name = "C2", show_name = F)
GSVA(gse_object = gse, gene_set = "C2", sub_set = "CGP", file_name = "C2_CGP", show_name = F)

