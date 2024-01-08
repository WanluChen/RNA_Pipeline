#################
## plot_PCA.R  ##
#################

library(ggplot2)
library(ggfortify)
library(DESeq2)
library(SummarizedExperiment)

gse <- readRDS("/Users/wanluchen/Documents/RNA_example_project/data/gse_object.rds")
gse$condition <- as.factor(gse$condition)

ddsTxi <- DESeqDataSet(gse, design = ~ condition)
dds <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 2,]
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) +
  labs(colour = "Condition")
# save plot in pdf
dir.create("/Users/wanluchen/Documents/RNA_example_project/results/PCA",
           showWarnings = FALSE, recursive = TRUE) ## change path here
pdf("/Users/wanluchen/Documents/RNA_example_project/results/PCA/PCA_group1.pdf",
    width = 5, height = 4) # change path here
p
dev.off()