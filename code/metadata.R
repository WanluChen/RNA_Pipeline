#################
## metadata.R  ##
#################

SraRunTable <- read.table(file = "/Users/wanluchen/Documents/RNA_example_project/data/SraRunTable.txt",
                          header = TRUE, sep = ",")
write.csv(SraRunTable, "/Users/wanluchen/Documents/RNA_example_project/data/SraRunTable.csv", row.names = FALSE)
