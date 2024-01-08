#####################
## process_data.R  ##
#####################

library(data.table)
library(dplyr)
library(GenomicRanges)
library(glue)
library(tximeta)
library(DBI)
library(biomaRt)
library(SummarizedExperiment)

# Metadata ---------------------------------------------------------------------
## change path here
samplesheet <- read_xlsx("/Users/wanluchen/Documents/RNA_example_project/data/samplesheet.xlsx")
samplesheet$files <- paste0("/Users/wanluchen/Documents/RNA_example_project/data/quants/", samplesheet$quants, "_quant/quant.sf")
write.csv(samplesheet, "/Users/wanluchen/Documents/RNA_example_project/data/sample_table.csv", row.names = F)

# Create SummarizedExperiment object -------------------------------------------
coldata <- read.csv("/Users/wanluchen/Documents/RNA_example_project/data/sample_table.csv")

se <- tximeta(coldata)
gse <- summarizeToGene(se)

# annotation -------------------------------------------------------------------
## choose your genome
mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl') ## human hg38
# mart <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl',
#                 host = "https://nov2020.archive.ensembl.org") ## mouse mm10
# mart <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl') ## mouse mm39

ens.str <- rownames(gse) <- substr(rownames(gse), 1, 15)
ensemble2gene <- getBM(attributes=c("ensembl_gene_id", "external_gene_name",
                                    "chromosome_name", "start_position", "end_position"),
                       filters = "ensembl_gene_id",
                       values = ens.str,
                       mart = mart)
rownames(ensemble2gene) <- ensemble2gene$ensembl_gene_id

rowData(gse)[rownames(ensemble2gene),"external_gene_name"] <- ensemble2gene$external_gene_name # duplicated with "gene_name"; if no "gene_name", use this
rowData(gse)[rownames(ensemble2gene),"chromosome_name"] <- ensemble2gene$chromosome_name
rowData(gse)[rownames(ensemble2gene),"start_position"] <- ensemble2gene$start_position
rowData(gse)[rownames(ensemble2gene),"end_position"] <- ensemble2gene$end_position

saveRDS(gse, "/Users/wanluchen/Documents/RNA_example_project/data/gse_object.rds")
