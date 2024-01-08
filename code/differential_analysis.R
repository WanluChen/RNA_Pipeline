##############################
## differential_analysis.R  ##
##############################

library(GenomicAlignments)
library(data.table)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(glue)


gse <- readRDS("/Users/wanluchen/Documents/RNA_example_project/data/gse_object.rds")

## if you want to change conditions, you can do it here
## for example
## gse$condition <- gsub("sg.*_","",gse$condition)

gse$condition <- as.factor(gse$condition)
ddsTxi <- DESeqDataSet(gse, design = ~ condition)
dds <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 2,]
dds <- DESeq(dds)


## run this section if you need GSEA later
if(1) {
  # generate expression file
  normalized_counts <- counts(dds, normalized = T)
  exp.dat <- data.frame(NAME = rowData(dds)$gene_name, DESCRIPTION = "NA")
  exp.dat <- cbind(exp.dat, normalized_counts)
  exp.dat <- exp.dat[order(exp.dat$NAME),]
  dir.create("/Users/wanluchen/Documents/RNA_example_project/results/gsea/input",
             showWarnings = FALSE, recursive = TRUE) ## change path here
  write.table(exp.dat, file = "/Users/wanluchen/Documents/RNA_example_project/results/gsea/input/expression.txt",
              sep = "\t", row.names = FALSE, quote = FALSE) ## change path here
  
  # n(samples) n(classes) 1
  # UserVisibleName_1 ... UserVisibleName_n(classes), (the order matches the unique names from the 3rd line)
  # Name_1 ... Name_n(samples), (classification for each sample)
  colnames(exp.dat)
  pheno <- rbind(paste(c(ncol(dds), length(unique(dds$condition)), 1), collapse = " "), ## change numbers here
                 paste0("# ", paste(unique(dds$condition), collapse = " ")), ## change names here
                 paste(dds$condition, collapse = " "))
  write.table(pheno, file = "/Users/wanluchen/Documents/RNA_example_project/results/gsea/input/ph.cls",
              sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE) ## change path here
}


# differential analysis
# Mutant vs Control
res <- results(dds, contrast=c("condition", "sg1_SMARCB1", "sgNT_Control"))
dir.create("/Users/wanluchen/Documents/RNA_example_project/results/differential_analysis/",
           showWarnings = FALSE, recursive = TRUE) ## change path here


# check the number of differential genes
## change threshold here
resOrdered <- res[order(res$padj),]
# significant: genes with adjusted p-value < 0.05 and log2 Fold Change > 0 or log2 Fold Change < 0
length(which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) > 0))
# up-regulated: genes with adjusted p-value < 0.05 and log2 Fold Change > 0
(n.up <- length(which(resOrdered$padj < 0.05 & resOrdered$log2FoldChange > 0)))
# down-regulated: genes with adjusted p-value < 0.05 and log2 Fold Change < 0
(n.down <- length(which(resOrdered$padj < 0.05 & resOrdered$log2FoldChange < 0)))


# Volcano plot
resOrdered <- data.frame(resOrdered)
de <- resOrdered[complete.cases(resOrdered), ] ## !!! NOTE: customize filtering here !!!
de$diffexpressed <- "Not Significant"
de$diffexpressed[de$log2FoldChange > 0 & de$padj < 0.1] <- glue("Up (n={n.up})")
de$diffexpressed[de$log2FoldChange < 0 & de$padj < 0.1] <- glue("Down (n={n.down})")
de$delabel <- NA
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  scale_color_manual(values=c("blue", "black", "red")) + 
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  labs(colour = "")
p
# save plot in pdf
pdf("/Users/wanluchen/Documents/RNA_example_project/results/differential_analysis/VolcanoPlot_pval0.05.pdf",
    width = 5, height = 4) ## change path here
p
dev.off()


# save results
resOrdered <- cbind(resOrdered, rowData(gse)[rownames(resOrdered),])
normalized_counts <- counts(dds, normalized = T) ## if you want normalized count, run this; otherwise skip
resOrdered <- cbind(resOrdered, normalized_counts[rownames(resOrdered),]) ## if you want normalized count, run this; otherwise skip
resOrdered <- apply(resOrdered,2,as.character)
write.csv(resOrdered, file = "/Users/wanluchen/Documents/RNA_example_project/results/differential_analysis/diffRes_sg1SMARCB1_sgNTControl.csv") ## change path here

