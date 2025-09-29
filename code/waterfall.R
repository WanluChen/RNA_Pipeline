library(ggplot2)
library(ggnewscale)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(glue)

get_ATAC_score <- function(gene_df, peak_df) {
  peak_df <- peak_df[!is.na(peak_df$SYMBOL) & !is.na(peak_df$atac_log2fc) & !is.na(peak_df$atac_padj), ]
  peak_df$gene <- peak_df$SYMBOL
  
  rna_atac <- merge(peak_df, gene_df, by = "gene")
  
  rna_atac$sign_match <- ifelse(sign(rna_atac$rna_log2fc) == sign(rna_atac$atac_log2fc), 1, -1)
  
  rna_atac$weighted_score <- rna_atac$sign_match * abs(rna_atac$atac_log2fc) * -log10(rna_atac$atac_padj)
  gene_score <- rna_atac %>%
    group_by(gene) %>%
    summarise(score_atac = sum(weighted_score), .groups = "drop")
  
  return(gene_score)
}

plot_and_save <- function(rna_df, atac_df, cond, cond2, filedir, suffix, tss_window = 2000, top_n = 50, filter_fun) {
  
  rna_sub <- rna_df[, c("gene_name", paste0("log2FoldChange_", cond), paste0("padj_", cond))]
  colnames(rna_sub) <- c("gene", "rna_log2fc", "rna_padj")
  rna_sub <- rna_sub[complete.cases(rna_sub) & rna_sub$gene != "", ]
  rna_sub <- rna_sub[rna_sub$rna_padj < 0.05, ]
  rna_sub$rna_padj[rna_sub$rna_padj == 0] <- 1e-300
  rna_sub$score_rna <- abs(rna_sub$rna_log2fc) * -log10(rna_sub$rna_padj)
  

  atac_sub <- atac_df[, c("chr", "start", "end", paste0("log2FoldChange_", cond), paste0("padj_", cond))]
  colnames(atac_sub) <- c("chr", "start", "end", "atac_log2fc", "atac_padj")
  atac_gr <- makeGRangesFromDataFrame(atac_sub, keep.extra.columns = TRUE)
  
  peak_annot <- annotatePeak(atac_gr,
                             TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                             annoDb = "org.Hs.eg.db",
                             tssRegion = c(-tss_window, tss_window),
                             verbose = FALSE)
  
  peak_df <- as.data.frame(peak_annot)
  peak_df <- peak_df[, c("seqnames", "start", "end", "SYMBOL", "distanceToTSS")]
  peak_df <- merge(peak_df, atac_sub, by.x = c("seqnames", "start", "end"), by.y = c("chr", "start", "end"))
  colnames(peak_df)[1] <- "chr"
  peak_df$atac_padj[peak_df$atac_padj == 0] <- 1e-300
  
  atac_score_df <- get_ATAC_score(gene_df = rna_sub, peak_df = peak_df)
  

  merged_score <- merge(rna_sub, atac_score_df, by = "gene")
  merged_score$combined_score <- merged_score$score_rna * merged_score$score_atac
  # merged_score$combined_score <- (merged_score$score_rna^0.8) * (merged_score$score_atac^0.2)
  
  merged_score <- merged_score[order(-merged_score$combined_score), ]
  merged_score <- merged_score[!duplicated(merged_score$gene), ]
  top_df <- filter_fun(merged_score)
  if (nrow(top_df) == 0) return(NULL)
  top_df <- top_df[1:min(top_n, nrow(top_df)), ]
  top_df <- top_df[order(top_df$rna_log2fc), ]
  top_df$gene <- factor(top_df$gene, levels = rev(top_df$gene))
  top_df$neglog10_padj <- -log10(top_df$rna_padj)
  top_df$clipped_rna_log2fc <- pmax(pmin(top_df$rna_log2fc, 5), -5)
  
  # ATAC bubble data
  bubble_data <- subset(peak_df, SYMBOL %in% top_df$gene)
  bubble_data$neglog10_padj <- -log10(bubble_data$atac_padj)
  bubble_data <- bubble_data[complete.cases(bubble_data[, c("atac_log2fc", "neglog10_padj")]), ]
  bubble_data$SYMBOL <- factor(bubble_data$SYMBOL, levels = levels(top_df$gene))
  bubble_data$scaled_atac_log2fc <- bubble_data$atac_log2fc
  
  # plot
  p <- ggplot() +
    geom_bar(data = top_df, aes(x = gene, y = clipped_rna_log2fc, fill = neglog10_padj), stat = "identity", alpha = 0.5) +
    scale_fill_gradient(low = "#F6E2FB", high = "#645A9F", name = "-log10(padj)") +
    
    ggnewscale::new_scale_fill() +
    
    geom_point(
      data = bubble_data,
      aes(x = SYMBOL, y = scaled_atac_log2fc,
          size = neglog10_padj, fill = atac_log2fc),
      shape = 21, color = "black", stroke = 0.3
    ) +
    scale_fill_gradient2(
      low = "#2E5A87", mid = "white", high = "#A90C38",
      midpoint = 0, name = "ATAC log2FC"
    ) +
    
    scale_size(range = c(2, 6)) +
    coord_flip() +
    labs(
      title = paste("Waterfall Plot:", cond, "vs", cond2, suffix),
      x = NULL, y = "log2 Fold Change", size = "-log10(padj)"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey80"),
      axis.ticks.x = element_blank()
    )
  
  # save
  outpath <- glue("results/waterfallplot/{filedir}/{cond}-{cond2}_{suffix}.pdf")
  ggsave(outpath, plot = p, width = 6, height = 8)
}

plot_and_save(rna_df, atac_df, cond, cond2 = "WT", filedir = "", suffix = "", top_n = 50, filter_fun = function(df) df)

