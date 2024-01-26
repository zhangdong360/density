library(readr)
library(tidyverse)
report_lib <- read_delim("data/report-lib.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
report_lib <- subset(report_lib,select = -c(FileName))
View(report_lib)
hist(report_lib$PrecursorMz)
plot(density(report_lib$PrecursorMz))

# 计算每个母离子对应的peptides
report_lib$PrecursorMz <- factor(report_lib$PrecursorMz)
result <- matrix(nrow = length(levels(report_lib$PrecursorMz)),
                 ncol = 3)
result <- as.data.frame(result)
colnames(result) <- c("precursorMz","peptide","protein")
# result <- data.frame(precursorMz = NA,
#                      peptide = NA,
#                      protein = NA
#                      )
for (i in 1:length(levels(report_lib$PrecursorMz))){
  report_one <- report_lib[report_lib$PrecursorMz == levels(report_lib$PrecursorMz)[i],]
  result$precursorMz[i] <- levels(report_lib$PrecursorMz)[i]
  result$peptide[i] <- length(unique(report_one$PeptideSequence))
  result$protein[i] <- length(unique(report_one$ProteinName))
  print(i)
}
# sum <- 0
# for (i in 1:33082){
#   if (result$peptide[i] == result$protein[1]){
#     sum <- sum + 1
#   }
# }
Mz_1 <- 1
Mz_2 <- 1
length(levels(report_lib$PrecursorMz))/160
z <- 25
for(i in 1:160){
  n <- i*z-z+1
  Mz_1[i] <- result$precursorMz[n]
  Mz_2[i] <- result$precursorMz[n+z]
}
report_lib$PrecursorMz <- as.numeric(report_lib$PrecursorMz)
Mz_1[1] <- 400
Mz_2[160] <- 1000
result_Mz <- data.frame(Mz_1 = as.numeric(Mz_1),
                        Mz_2 = as.numeric(Mz_2))
result_Mz$windows_size <- result_Mz$Mz_2-result_Mz$Mz_1
result_Mz$`m/z range` <- paste0(result_Mz$Mz_1,"-",result_Mz$Mz_2)
write.csv(result_Mz,file = 'result_Mz.csv')
