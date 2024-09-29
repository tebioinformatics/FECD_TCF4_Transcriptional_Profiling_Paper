### Heatmap Analysis
# This script performs heatmap analysis using TPM expression data and clusters both genes and samples.

# Clear all previous data from the environment
rm(list=ls())

# Set the working directory (replace with your actual path)
working_directory <- "path/to/your/working/directory"
setwd(working_directory)

# Load the expression data (e.g., TPM values) from the file
expression_file <- "DataMatrix_ImportSampleList_Case_vs_Control.txt"  # Replace with your actual file
expression_data <- read.table(expression_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Check dimensions and preview the data
dim(expression_data)
head(expression_data)

# Extract abundance data (first 10 columns, in this example)
abundance_data <- expression_data[, 1:10]  # Modify this based on the actual number of samples
head(abundance_data)

# Load the sample metadata file (containing sample names and conditions)
sample_list_file <- "Sample_List_YYYYMMDD.txt"  # Replace with your actual file
sample_metadata <- read.table(sample_list_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set the column names of the abundance data to the sample names from the sample metadata
colnames(abundance_data) <- sample_metadata$sample
head(abundance_data)

# Filter genes with TPM >= 1 in all samples
library(genefilter)  # For filtering gene expression data
filter_condition <- kOverA(10, A = 1)  # Filter condition: genes with TPM >= 1 in 10 or more samples
filter_function <- filterfun(filter_condition)  # Create filter function
filtered_genes <- genefilter(abundance_data, filter_function)  # Apply the filter
filtered_abundance <- abundance_data[filtered_genes, ]  # Subset the data based on the filter
dim(filtered_abundance)

# Apply log10 transformation to the filtered data
log_transformed_data <- log10(filtered_abundance + 1)
dim(log_transformed_data)

### Create the heatmap
heatmap_matrix <- as.matrix(log_transformed_data)

# Calculate distances for clustering
gene_distance <- dist(heatmap_matrix, method = "manhattan")  # Distance between genes
sample_distance <- dist(t(heatmap_matrix), method = "manhattan")  # Distance between samples
# You can change the distance method: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"

# Perform hierarchical clustering
gene_clustering <- hclust(gene_distance, method = "ward.D2")
sample_clustering <- hclust(sample_distance, method = "ward.D2")
# You can change the clustering method: "ward.D", "ward.D2", "single", "complete", "average", etc.

# Load required libraries for plotting the heatmap
library(gplots)  # For creating heatmaps
library(RColorBrewer)  # For color palettes

# Plot the heatmap with row and column clustering
heatmap.2(heatmap_matrix,
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),  # Color palette (Red to Blue)
          scale = "row",  # Scale data by row (z-score normalization across genes)
          key = TRUE,  # Show color key
          keysize = 1.0,  # Size of the color key
          density.info = "none",  # No density info
          trace = "none",  # No trace lines
          cexCol = 1,  # Font size for column labels (samples)
          cexRow = 0.3,  # Font size for row labels (genes)
          margins = c(3, 10),  # Margins for the heatmap plot
          Rowv = as.dendrogram(gene_clustering),  # Apply hierarchical clustering to rows (genes)
          Colv = as.dendrogram(sample_clustering)  # Apply hierarchical clustering to columns (samples)
)

# Optionally save the heatmap plot (you can adjust the file path and format as needed)
heatmap_file <- "Heatmap_Case_vs_Control.png"  # Replace with your desired output file name
dev.copy(png, file = heatmap_file, width = 800, height = 600)
dev.off()
