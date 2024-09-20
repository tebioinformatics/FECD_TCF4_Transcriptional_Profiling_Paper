# 2024/08/29
# Differential gene expression testing from expression data using DESeq2
# Remove all previous data


#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install()

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("DESeq2")

#source ("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#library("DESeq2")


#source ("http://bioconductor.org/biocLite.R")

#BiocManager::install("tximport")

rm(list=ls())
##############################################################################################

path <- "/Users/shara/Library/CloudStorage/OneDrive-同志社大学/遺伝子解析班データファイル/RNA-Seqプロジェクト/DESeq2/DU2022_QC4_v104/gene_RE-_vs_Control_20240829AT"
setwd(path)

library(DESeq2)
library(tximport)

# Load the list containing sample information
s2c <- read.table("Import_Sample_List_RE-_vs_Control_20240829AT.txt", header=T, sep="\t", stringsAsFactors=F)
files <- s2c$path
names(files) <- s2c$sample

# Load expression data (analysis at gene-level)
txi <- tximport(files, type="rsem", txIn=F, txOut=F)
dim(as.data.frame(txi))
write.table(txi, file="DataMatrix_ImportSampleList_RE-_vs_Control_20240829AT.txt", sep="\t", quote=F)

# Replace length=0 with 1
txi$length[txi$length==0] <- 1
sampleTable <- data.frame(condition=s2c$group,RIN=s2c$RIN)
rownames(sampleTable) <- colnames(txi$counts)
# Explicitly specify factor type
sampleTable$condition <- as.factor(sampleTable$condition)　# Explicitly specify factor type


dds <- DESeqDataSetFromTximport(txi, sampleTable,
                                design = ~ RIN + condition)



# Wald test (2-group comparison)
dds_wt <- DESeq(dds)
res_wt <- results(dds_wt, contrast=c("condition", "nex","control")) 　# Specify the order of subtraction for log2FoldChange
dim(res_wt)

summary(res_wt)

res_wt_naomit <- na.omit(res_wt) # Remove missing values
dim(res_wt_naomit)

# Write out data
write.table(res_wt_naomit, file="DESeq2_Result_of_WaldTest_RE-_vs_Control_20240829AT.txt", sep="\t", quote=F)
