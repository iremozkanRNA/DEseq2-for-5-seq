#Irem 5seq DEseq2analysis prep
rm(list=ls())
library(dplyr)

files <- list.files(path ="~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/5seqBeds/5seqBedsOnly/Ksg_NDC/Processed/" , pattern = "*.bed",full.names = T)
temp <- read.delim(files[[1]], header = FALSE)
temp1 <- read.delim(files[[2]], header = FALSE)
temp3 <- read.delim(files[[3]], header = F)
outputfilename <- "~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/5seqBeds/5seqBedsOnly/CountMatrix/316-Processed_KsgNDCt60_ReadLenght_CountMatrix.csv"
analyzed <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/5seq/TSSPredIremOutput/KsgData_nofRank.csv") %>% filter(Genome == "NDCt60")

temp <- temp[,-c(4,5)]

tempMinus <- temp[temp$V6 %in% "-",]
tempPlus <- temp[temp$V6 %in% "+",]

tempMinusCounts <- as.data.frame(table(tempMinus$V2))

tempPlusCounts <- as.data.frame(table(tempPlus$V3))


tempMinus <- distinct(tempMinus, V2, .keep_all = TRUE)
tempMinus <- arrange(tempMinus, V2)
tempMinus$coverage <- tempMinusCounts$Freq[match(tempMinus$V2, tempMinusCounts$Var1)]

tempPlus <- distinct(tempPlus, V3, .keep_all = TRUE)
tempPlus <- arrange(tempPlus, V3)
tempPlus$coverage <- tempPlusCounts$Freq[match(tempPlus$V3, tempPlusCounts$Var1)]


tempMinus <- tempMinus[,c(1,2,5,4)]
colnames(tempMinus) <- c("chrom","readPosition","coverage","strand")

tempPlus <- tempPlus[,c(1,3,5,4)]
colnames(tempPlus) <- c("chrom","readPosition","coverage","strand")

tempOut <- rbind(tempMinus, tempPlus)
tempOut <- tempOut %>% arrange(readPosition)

tempo <- temp %>% distinct(V3, .keep_all = T)
tempOut$Readstart <- ""
tempOut$Readend <- ""

matches <- match(tempOut$readPosition, tempo$V3)
tempOut$Readstart <- tempo$V2[matches]
tempOut$Readend <- tempOut$readPosition

tempOutpos <- tempOut %>% filter(strand=="+")
tempoutneg <- tempOut %>% filter(strand=="-")


analyzedNeg <- analyzed %>% filter(SuperStrand == "-")
analyzedPos <- analyzed %>% filter(SuperStrand == "+")

countMatrixPos <- data.frame(
  chrom = "NC_003028",
  SuperPos = analyzedPos$SuperPos,
  HighestPeakCoverage_a = ""
)
# Iterate over each SuperPos
for (i in seq_along(countMatrixPos$SuperPos)) {
  current_pos <- countMatrixPos$SuperPos[i]
  previous_pos <- ifelse(i == 1, -Inf, countMatrixPos$SuperPos[i - 1])
  
  # Filter temput for positions within the current range
  filtered_reads <- tempOutpos %>%
    filter(readPosition > (current_pos -54) & (readPosition <= current_pos+54))
  
  # Calculate total coverage for these reads
  total_coverage <- sum(filtered_reads$coverage)
  
  # Store results
  countMatrixPos$HighestPeakCoverage_a[i] <- total_coverage
}

countMatrixNeg <- data.frame(
  chrom = "NC_003028",
  SuperPos = analyzedNeg$SuperPos,
  HighestPeakCoverage_a = ""
)

# Iterate over each SuperPos
for (i in seq_along(countMatrixNeg$SuperPos)) {
  current_pos <- countMatrixNeg$SuperPos[i]
  previous_pos <- ifelse(i == 1, -Inf, countMatrixNeg$SuperPos[i - 1])
  
  # Filter temput for positions within the current range
  filtered_reads <- tempOutpos %>%
    filter(readPosition > (current_pos -54) & (readPosition <= current_pos+54))
  
  # Calculate total coverage for these reads
  total_coverage <- sum(filtered_reads$coverage)
  
  # Store results
  countMatrixNeg$HighestPeakCoverage_a[i] <- total_coverage
}

CountMatrix <- rbind(countMatrixPos, countMatrixNeg)
CountMatrix <- CountMatrix %>% arrange(SuperPos)
CountMatrix <- CountMatrix[!duplicated(CountMatrix$SuperPos),]
rownames(CountMatrix) <- CountMatrix$SuperPos
CountMatrix <- CountMatrix[, -1] # Remove region column
CountMatrixMeta <- CountMatrix

######################## Process file b

temp1 <- temp1[,-c(4,5)]

temp1Minus <- temp1[temp1$V6 %in% "-",]
temp1Plus <- temp1[temp1$V6 %in% "+",]

temp1MinusCounts <- as.data.frame(table(temp1Minus$V2))

temp1PlusCounts <- as.data.frame(table(temp1Plus$V3))


temp1Minus <- distinct(temp1Minus, V2, .keep_all = TRUE)
temp1Minus <- arrange(temp1Minus, V2)
temp1Minus$coverage <- temp1MinusCounts$Freq[match(temp1Minus$V2, temp1MinusCounts$Var1)]

temp1Plus <- distinct(temp1Plus, V3, .keep_all = TRUE)
temp1Plus <- arrange(temp1Plus, V3)
temp1Plus$coverage <- temp1PlusCounts$Freq[match(temp1Plus$V3, temp1PlusCounts$Var1)]


temp1Minus <- temp1Minus[,c(1,2,5,4)]
colnames(temp1Minus) <- c("chrom","readPosition","coverage","strand")

temp1Plus <- temp1Plus[,c(1,3,5,4)]
colnames(temp1Plus) <- c("chrom","readPosition","coverage","strand")

temp1Out <- rbind(temp1Minus, temp1Plus)
temp1Out <- temp1Out %>% arrange(readPosition)

temp1o <- temp1 %>% distinct(V3, .keep_all = T)
temp1Out$Readstart <- ""
temp1Out$Readend <- ""

matches <- match(temp1Out$readPosition, temp1o$V3)
temp1Out$Readstart <- temp1o$V2[matches]
temp1Out$Readend <- temp1Out$readPosition

temp1Outpos <- temp1Out %>% filter(strand=="+")
temp1outneg <- temp1Out %>% filter(strand=="-")

countMatrixPos <- data.frame(
  chrom = "NC_003028",
  SuperPos = analyzedPos$SuperPos,
  HighestPeakCoverage = ""
)
# Iterate over each SuperPos
for (i in seq_along(countMatrixPos$SuperPos)) {
  current_pos <- countMatrixPos$SuperPos[i]
  previous_pos <- ifelse(i == 1, -Inf, countMatrixPos$SuperPos[i - 1])
  
  # Filter temp1ut for positions within the current range
  filtered_reads <- temp1Outpos %>%
    filter(readPosition > (current_pos -54) & (readPosition <= current_pos+54))
  
  # Calculate total coverage for these reads
  total_coverage <- sum(filtered_reads$coverage)
  
  # Store results
  countMatrixPos$HighestPeakCoverage[i] <- total_coverage
}

countMatrixNeg <- data.frame(
  chrom = "NC_003028",
  SuperPos = analyzedNeg$SuperPos,
  HighestPeakCoverage = ""
)

# Iterate over each SuperPos
for (i in seq_along(countMatrixNeg$SuperPos)) {
  current_pos <- countMatrixNeg$SuperPos[i]
  previous_pos <- ifelse(i == 1, -Inf, countMatrixNeg$SuperPos[i - 1])
  
  # Filter temp1ut for positions within the current range
  filtered_reads <- temp1Outpos %>%
    filter(readPosition > (current_pos -54) & (readPosition <= current_pos+54))
  
  # Calculate total coverage for these reads
  total_coverage <- sum(filtered_reads$coverage)
  
  # Store results
  countMatrixNeg$HighestPeakCoverage[i] <- total_coverage
}

CountMatrix <- rbind(countMatrixPos, countMatrixNeg)
CountMatrix <- CountMatrix %>% arrange(SuperPos)
CountMatrix <- CountMatrix[!duplicated(CountMatrix$SuperPos),]
rownames(CountMatrix) <- CountMatrix$SuperPos
CountMatrix <- CountMatrix[, -1] # Remove region column

CountMatrixMeta$HighestPeakCoverage_b <- CountMatrix$HighestPeakCoverage



######################## Process file c

temp2 <- temp2[,-c(4,5)]

temp2Minus <- temp2[temp2$V6 %in% "-",]
temp2Plus <- temp2[temp2$V6 %in% "+",]

temp2MinusCounts <- as.data.frame(table(temp2Minus$V2))

temp2PlusCounts <- as.data.frame(table(temp2Plus$V3))


temp2Minus <- distinct(temp2Minus, V2, .keep_all = TRUE)
temp2Minus <- arrange(temp2Minus, V2)
temp2Minus$coverage <- temp2MinusCounts$Freq[match(temp2Minus$V2, temp2MinusCounts$Var1)]

temp2Plus <- distinct(temp2Plus, V3, .keep_all = TRUE)
temp2Plus <- arrange(temp2Plus, V3)
temp2Plus$coverage <- temp2PlusCounts$Freq[match(temp2Plus$V3, temp2PlusCounts$Var1)]


temp2Minus <- temp2Minus[,c(1,2,5,4)]
colnames(temp2Minus) <- c("chrom","readPosition","coverage","strand")

temp2Plus <- temp2Plus[,c(1,3,5,4)]
colnames(temp2Plus) <- c("chrom","readPosition","coverage","strand")

temp2Out <- rbind(temp2Minus, temp2Plus)
temp2Out <- temp2Out %>% arrange(readPosition)

temp2o <- temp2 %>% distinct(V3, .keep_all = T)
temp2Out$Readstart <- ""
temp2Out$Readend <- ""

matches <- match(temp2Out$readPosition, temp2o$V3)
temp2Out$Readstart <- temp2o$V2[matches]
temp2Out$Readend <- temp2Out$readPosition

temp2Outpos <- temp2Out %>% filter(strand=="+")
temp2outneg <- temp2Out %>% filter(strand=="-")

countMatrixPos <- data.frame(
  chrom = "NC_003028",
  SuperPos = analyzedPos$SuperPos,
  HighestPeakCoverage = ""
)
# Iterate over each SuperPos
for (i in seq_along(countMatrixPos$SuperPos)) {
  current_pos <- countMatrixPos$SuperPos[i]
  previous_pos <- ifelse(i == 1, -Inf, countMatrixPos$SuperPos[i - 1])
  
  # Filter temp2ut for positions within the current range
  filtered_reads <- temp2Outpos %>%
    filter(readPosition > (current_pos -54) & (readPosition <= current_pos+54))
  
  # Calculate total coverage for these reads
  total_coverage <- sum(filtered_reads$coverage)
  
  # Store results
  countMatrixPos$HighestPeakCoverage[i] <- total_coverage
}

countMatrixNeg <- data.frame(
  chrom = "NC_003028",
  SuperPos = analyzedNeg$SuperPos,
  HighestPeakCoverage = ""
)

# Iterate over each SuperPos
for (i in seq_along(countMatrixNeg$SuperPos)) {
  current_pos <- countMatrixNeg$SuperPos[i]
  previous_pos <- ifelse(i == 1, -Inf, countMatrixNeg$SuperPos[i - 1])
  
  # Filter temp2ut for positions within the current range
  filtered_reads <- temp2Outpos %>%
    filter(readPosition > (current_pos -54) & (readPosition <= current_pos+54))
  
  # Calculate total coverage for these reads
  total_coverage <- sum(filtered_reads$coverage)
  
  # Store results
  countMatrixNeg$HighestPeakCoverage[i] <- total_coverage
}

CountMatrix <- rbind(countMatrixPos, countMatrixNeg)
CountMatrix <- CountMatrix %>% arrange(SuperPos)
CountMatrix <- CountMatrix[!duplicated(CountMatrix$SuperPos),]
rownames(CountMatrix) <- CountMatrix$SuperPos
CountMatrix <- CountMatrix[, -1] # Remove region column

CountMatrixMeta$HighestPeakCoverage_c <- CountMatrix$HighestPeakCoverage

write.csv(CountMatrixMeta, file=outputfilename, col.names = T, row.names = F)
