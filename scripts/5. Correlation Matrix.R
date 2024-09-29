### Correlation Matrix Analysis
# This script performs correlation matrix analysis using gene expression data (TPM) and visualizes the result.

# Clear all previous data
rm(list=ls())

# Set the working directory (replace with your actual path)
working_directory <- "path/to/your/working/directory"
setwd(working_directory)

# Load expression data (e.g., TPM values)
tpm_file <- "DataMatrix_ImportSampleList_Case_vs_Control.txt"  # Replace with your actual file
tpm_data <- read.table(tpm_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract abundance data (e.g., TPM values for 10 samples)
abundance_data <- tpm_data[, 1:10]
colnames(abundance_data)
head(abundance_data)

# Load the sample metadata (containing sample names and conditions)
sample_metadata_file <- "Sample_List_YYYYMMDD.txt"  # Replace with your actual file
sample_metadata <- read.table(sample_metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set the column names of abundance data to sample names from the sample metadata
colnames(abundance_data) <- sample_metadata$sample
head(abundance_data)
dim(abundance_data)

# Filter genes with TPM >= 1 in all samples
library(genefilter)  # For filtering gene expression data
filter_condition <- kOverA(10, A = 1)  # Filter condition: genes with TPM >= 1 in 10 or more samples
filter_function <- filterfun(filter_condition)  # Create a filtering function
filtered_genes <- genefilter(abundance_data, filter_function)  # Apply the filter
filtered_tpm <- abundance_data[filtered_genes, ]  # Subset the data based on the filter
dim(filtered_tpm)

# Apply log10 transformation to the filtered TPM data
log_transformed_data <- log10(filtered_tpm + 1)
dim(log_transformed_data)

### Calculate Correlation Matrix
# Compute Spearman correlation coefficients between samples
correlation_matrix <- cor(log_transformed_data, method = "spearman")

# Round the correlation values to two decimal places for readability
rounded_correlation_matrix <- round(correlation_matrix, digits = 2)

# Load the required library for correlation matrix visualization
library(corrplot)
packageVersion("corrplot")  # Check the package version

# Define the color palette for the correlation matrix
color_palette <- colorRampPalette(c("blue", "cyan", "white", "yellow", "red"))

### Plot the Correlation Matrix
corrplot(rounded_correlation_matrix,
         method = "color",  # Use color shading to represent the correlation values
         shade.col = NA,  # No shading color
         tl.col = "black",  # Color for text labels (sample names)
         number.cex = 1.5,  # Font size for correlation coefficients
         addCoef.col = "black",  # Color for correlation coefficient numbers
         tl.cex = 1.2,  # Font size for text labels
         cl.cex = 1.2,  # Font size for color legend
         order = "hclust",  # Order samples by hierarchical clustering
         col = color_palette(200),  # Color palette for the heatmap
         tl.srt = 45,  # Rotate sample names 45 degrees for better readability
         col.lim = c(-1, 1)  # Set color scale limits to [-1, 1] (full correlation range)
)
