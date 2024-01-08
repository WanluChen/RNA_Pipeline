###############
## get_QC.R  ##
###############

library(glue)
library(data.table)
library(magrittr)
library(rjson)

filename <- dir(path = "/Users/wanluchen/Documents/RNA_example_project/data/quants", full.names = TRUE) ## change path here

mapping_rate_all <- list()
for(i in filename){
  
  data <- fromJSON(file = glue("{i}/aux_info/meta_info.json"))
  
  mapping_rate <- data$percent_mapped
  total_read <- data$num_processed
  map_read <- data$num_mapped
  sample_name <- sub(".*/","",i) %>% sub("_quant.*","",.)
  
  mapping_rate_all[[i]] <- data.frame(quants = sample_name, mapping_rate = mapping_rate,
                                      sequence_read = total_read, map_read = map_read)
}
mapping_rate_combine <- Reduce(rbind,mapping_rate_all)
align_info <- data.frame(mapping_rate_combine)

samplesheet <- read_excel("/Users/wanluchen/Documents/RNA_example_project/data/samplesheet.xlsx")
df <- merge(samplesheet, align_info, by = "quants")

dir.create("/Users/wanluchen/Documents/RNA_example_project/results/QC", showWarnings = FALSE, recursive = TRUE) ## change path here
write.csv(df, file="/Users/wanluchen/Documents/RNA_example_project/results/QC/mapping_rate_rna.csv", row.names = FALSE) ## change path here
