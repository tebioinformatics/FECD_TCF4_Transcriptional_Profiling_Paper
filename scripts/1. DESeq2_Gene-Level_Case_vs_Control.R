# Differential gene expression testing from expression data using DESeq2
# This script performs differential expression analysis on RNA-Seq data using the DESeq2 package

# Clear all objects from the environment to ensure a fresh start
rm(list=ls())

# Set the working directory (please adjust the path as necessary)
working_directory <- "path/to/your/working/directory"
setwd(working_directory)

# Load required libraries
# If DESeq2 or tximport are not installed, uncomment and run the following lines:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")

library(DESeq2)   # For differential expression analysis
library(tximport) # For importing transcript-level estimates from different quantification tools


# Load the list containing sample information
# The file should contain columns for sample paths, sample names, conditions, and any other metadata such as RIN
sample_info_file <- "Sample_List_YYYYMMDD.txt"
s2c <- read.table(sample_info_file, header=T, sep="\t", stringsAsFactors=F)
files <- s2c$path
names(files) <- s2c$sample

# Load expression data (analysis at gene-level)
txi <- tximport(files, type="rsem", txIn=F, txOut=F)
dim(as.data.frame(txi))
output_matrix_file <- "DataMatrix_YYYYMMDD.txt"
write.table(txi, file=output_matrix_file, sep="\t", quote=F)

# Replace length=0 with 1
txi$length[txi$length==0] <- 1
sampleTable <- data.frame(condition=s2c$group,RIN=s2c$RIN)
rownames(sampleTable) <- colnames(txi$counts)

# Explicitly specify factor type
sampleTable$condition <- as.factor(sampleTable$condition)

# Create a DESeq2 dataset from the tximport object and sample table
# The design formula accounts for the effects of RNA quality (RIN) and experimental condition
dds <- DESeqDataSetFromTximport(txi, sampleTable,
                                design = ~ RIN + condition)

# Perform differential expression analysis using the Wald test (2-group comparison)
# Here, we compare "case" vs. "control" (ctr). Adjust the contrast as necessary for your study.
dds_wt <- DESeq(dds)
res_wt <- results(dds_wt, contrast=c("condition", "case","ctr")) # Specify the order of subtraction for log2FoldChange
dim(res_wt)

# Summarize the results of the Wald test
summary(res_wt)

# Remove rows with missing values (NA) from the results
res_wt_naomit <- na.omit(res_wt)
dim(res_wt_naomit)

# Save the filtered results to a file for further analysis or reporting
output_file <- "DESeq2_Result_of_WaldTest_Case_vs_Control.txt"
write.table(res_wt_naomit, file=output_file, sep="\t", quote=F)
