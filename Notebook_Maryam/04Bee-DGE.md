# Differential gene expression with DeSeq2


This code should be run in R:

```
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")

dat<-read.table("At_count.txt",header = T,quote = "",row.names = 1)

# Convert to matrix
dat <- as.matrix(dat)
