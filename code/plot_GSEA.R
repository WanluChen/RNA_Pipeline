library(ggplot2)
library(ggtext)
library(dplyr)
library(readr)
library(glue)
library(stringr)


barplotGSEA <- function(dir1, dir2 = NULL, cond1, cond2, geneset) {
  if (is.null(dir2)) dir2 <- paste(cond1, cond2, geneset, sep = "_")
  
  # read
  dirs <- list.dirs(path = glue("results/gsea/{dir1}"), full.names = TRUE, recursive = FALSE)
  filtered_dirs <- dirs[grep(paste0("^", dir2), basename(dirs))]
  files <- list.files(path = filtered_dirs, pattern = "^gsea_report_for.*\\.tsv$", full.names = FALSE)
  
  gsea1 <- read_tsv(file.path(filtered_dirs, files[grep(paste0("^gsea_report_for_", cond1, "_"), files)]), show_col_types = FALSE)
  gsea2 <- read_tsv(file.path(filtered_dirs, files[grep(paste0("^gsea_report_for_", cond2, "_"), files)]), show_col_types = FALSE)
  gsea1$Source <- cond1
  gsea2$Source <- cond2
  gsea_combined <- bind_rows(gsea1, gsea2)
  
  # FDR < 0.25 & top 20
  top_gsea <- gsea_combined %>%
    rename(FDR = `FDR q-val`) %>%
    filter(FDR < 0.25) %>%
    group_by(Source) %>%
    arrange(FDR) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  # plot
  top_gsea$NAME_wrapped <- str_wrap(gsub("_", " ", top_gsea$NAME), width = 50)
  p <- ggplot(top_gsea, aes(x = NES, y = reorder(NAME_wrapped, NES), fill = FDR)) +
    geom_col() +
    scale_fill_gradient(low = "#4E7CA1", high = "#C3D7DF") +
    labs(
      title = "GSEA Bar Plot (Top 20, FDR < 0.25)",
      subtitle = glue("**{cond1}** vs **{cond2}** on GeneSet **{toupper(geneset)}**"),
      x = "Normalized Enrichment Score (NES)",
      y = "Gene Set",
      fill = "FDR"
    ) +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title = element_markdown(hjust = 0),
      plot.subtitle = element_markdown(hjust = 0),
      axis.text.y = element_text(size = 15)
    )
  
  # save
  out_dir <- glue("results/gsea/{dir1}/top20_fdr0.25")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(p, filename = glue("{out_dir}/{dir2}_barplot.pdf"), width = 12, height = 15)
}

comparisons <- data.frame(
  cond1 = c("PRCC_TFE3", "PRCCdPRW_TFE3", "TFE3"),
  cond2 = c("TFE3", "PRCC_TFE3", "KO")
)

for (i in 1:nrow(comparisons)) {
  for (geneset in c("c2", "c5", "h")) {
    barplotGSEA(dir1 = "v1_h+c2+c5", cond1 = comparisons$cond1[i], cond2 = comparisons$cond2[i], geneset = geneset)
  }
}


