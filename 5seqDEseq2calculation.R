rm(list=ls())
library(data.table)
library(dplyr)
setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/5seqBeds/5seqBedsOnly/CountMatrix/")
outputfile <- "~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/5seqBeds/5seqDEseq2/Only54eachend/DEseq2_Chl3QvsNDCresults.csv"
# Load counts from multiple samples
countsMeta <- read.csv("316-Processed_ChlNDCt60_ReadLenght_CountMatrix.csv")
#counts2 <- read.csv("316-UnProcessed_KsgNDCt60_ReadLenght_CountMatrix.csv")
#countsMeta <- counts
#countsMeta[,c(5:7)] <- counts2[,c(2:4)]
colnames(countsMeta) <- c("SuperPos","ctr_Prep_a", "ctr_Prep_b", "ctr_Prep_c")


# Enter Your comparison Files
countsMeta2 <- read.csv("316-Processed_Chl3qt60_ReadLenght_CountMatrix.csv")
#counts4 <- read.csv("316-UnProcessed_Ksg3qt60_ReadLenght_CountMatrix.csv")
#countsMeta2 <- counts3
#countsMeta2[,c(5:7)] <- counts4[,c(2:4)]
colnames(countsMeta2) <- c("SuperPos","Prep_a", "Prep_b", "Prep_c")

countsMeta3 <- full_join(countsMeta, countsMeta2, by = "SuperPos")
countsMeta3[is.na(countsMeta3)] <- 0
rownames(countsMeta3) <- countsMeta3$SuperPos
countsMeta3 <- countsMeta3[,-c(1)]

sampleTable <- data.frame(
  sampleName = c("ctr_Prep_a", "ctr_Prep_b", "ctr_Prep_c", "Prep_a", "Prep_b", "Prep_c"),
  condition = factor(c("ctr","ctr", "ctr","trt","trt","trt"))
)


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = countsMeta3,
                              colData = sampleTable,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
rest <- as.data.frame(res)
# Save results to file
write.csv(as.data.frame(res), outputfile, col.names = T, row.names = T)
