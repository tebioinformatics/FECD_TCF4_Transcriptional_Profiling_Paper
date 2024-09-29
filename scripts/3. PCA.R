### Principal Component Analysis (PCA)
# This script performs a PCA analysis using TPM (transcripts per million) expression data and visualizes the result.

# Clear all previous data
rm(list=ls())

# Set the working directory (replace with your actual path)
working_directory <- "path/to/your/working/directory"
setwd(working_directory)

# Load required libraries
library(genefilter)  # For filtering genes based on expression levels
library(ggplot2)     # For creating plots

# Load expression data (TPM: Transcripts Per Million)
tpm_file <- "DataMatrix_ImportSampleList_Case_vs_Control.txt" #(For exaple, 4 RE- and 6 RE+ in this case)
tpm_data <- read.table(tpm_file, sep = "\t", header = TRUE)

# Extract TPM data (first 10 columns as gene expression levels for the samples)
tpm_matrix <- tpm_data[, 1:10]
head(tpm_matrix)
dim(tpm_matrix)

# Filter genes with TPM >= 1 in all samples
filter_condition <- kOverA(10, A = 1)  # Filter for genes with TPM >= 1 in 10 or more samples
filter_function <- filterfun(filter_condition)  # Create a filtering function
filtered_genes <- genefilter(tpm_matrix, filter_function)  # Apply filter to extract genes that meet the condition
filtered_tpm <- tpm_matrix[filtered_genes, ]  # Filtered TPM data

# Check dimensions after filtering
dim(filtered_tpm)
head(filtered_tpm)

# Load the sample information list (metadata for the samples)
sample_list_file <- "Sample_List_YYYYMMDD.txt"
sample_list <- read.table(sample_list_file, sep = "\t", header = TRUE)

# Set column names for TPM data to sample names from the sample list
colnames(filtered_tpm) <- sample_list$sample
head(filtered_tpm)

### PCA Calculation
# Apply log10 transformation to TPM data
log_tpm <- log10(filtered_tpm + 1)

# Perform PCA (Principal Component Analysis)
pca_result <- prcomp(t(log_tpm), scale = FALSE)  # Transpose the matrix for PCA; genes are in rows, samples in columns

# Check the proportion of variance explained by the first two principal components
summary(pca_result)$importance

# Extract the first two principal components (PC1 and PC2)
pc1_values <- pca_result$x[, 1]
pc2_values <- pca_result$x[, 2]

# Create a data frame for PCA plot
pca_data <- data.frame(PC1 = pc1_values, PC2 = pc2_values)
sample_groups <- c(rep("RE-", 4), rep("RE+", 6))  # Adjust based on your actual grouping (4 RE-, 6 RE+ in this case)
sample_names <- colnames(log_tpm)
pca_data <- cbind(Sample = sample_names, pca_data, Group = sample_groups)

# Calculate percentage variance for PC1 and PC2
variance_pc1 <- signif(summary(pca_result)$importance[2, 1] * 100, 3)
variance_pc2 <- signif(summary(pca_result)$importance[2, 2] * 100, 3)

# Labels for PC axes with percentage variance
pc1_label <- paste("PC1 (", variance_pc1, "%)", sep = "")
pc2_label <- paste("PC2 (", variance_pc2, "%)", sep = "")

### PCA Plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, group = Group, fill = Group)) +
  theme_linedraw() +  # Simple theme with lines
  labs(x = pc1_label, y = pc2_label) +  # Axis labels
  geom_text(aes(label = Sample), vjust = -1, size = 5) +  # Sample names as labels on the plot
  geom_point(aes(shape = Group), size = 7, stroke = 0.8) +  # Points for each sample, differentiated by shape
  scale_shape_manual(values = c("RE-" = 23, "RE+" = 21)) +  # Shape for groups
  scale_fill_manual(values = c("#67A9CF", "#EF8A62")) +  # Colors for groups
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1.0),  # Black border around the plot
    panel.grid = element_line(size = 1.0, colour = "black", linetype = "dashed"),  # Grid lines
    aspect.ratio = 1.0,  # Aspect ratio of the plot
    text = element_text(size = 17)  # Font size for labels
  ) +
  xlim(-20, 20) +  # Adjust x-axis limits
  ylim(-20, 20)  # Adjust y-axis limits

# Display the plot
print(pca_plot)

# Optionally save the plot as a PNG file
ggsave(filename = "PCA_Plot_RE-_vs_RE+.png", 
       plot = pca_plot, 
       device = "png", 
       width = 6, height = 6, 
       units = "in", 
       dpi = 300)
