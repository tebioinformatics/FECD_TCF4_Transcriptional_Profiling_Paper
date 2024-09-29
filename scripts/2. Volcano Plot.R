### Volcano Plot for Differential Expression Analysis
# This script creates a volcano plot based on DESeq2 differential expression results.

# Clear environment and close open plots (if any)
rm(list=ls()) 
#dev.off() # Uncomment if you need to close any open graphics window

### Set working directory and read data
working_directory <- "path/to/your/working/directory/" # Replace with your directory path
setwd(working_directory)

# Load the DESeq2 results
deseq2_result_file <- "DESeq2_Result_of_WaldTest_Case_vs_Control.txt"
deseq2_data <- read.table(deseq2_result_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract relevant columns (padj and log2FoldChange) and add to a new data frame
volcano_data <- cbind(padj = deseq2_data$padj, log2FoldChange = deseq2_data$log2FoldChange)
rownames(volcano_data) <- rownames(deseq2_data)
head(volcano_data)

# Convert p-adjusted values to -log10(p.adjust)
volcano_data[,1] <- -log10(volcano_data[,1])  # Convert p.adjust for plotting
head(volcano_data)

# Load additional expression level data (TPM)
tpm_directory <- "path/to/your/working/directory/" # Replace with your directory path
setwd(tpm_directory)

tpm_file <- "DataMatrix_ImportSampleList_Case_vs_Control.txt"  # Input file for TPM
tpm_data <- read.table(tpm_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter genes based on average expression level (TPM >= 1)
filtered_tpm_data <- tpm_data[rowMeans(tpm_data[, 1:17]) >= 1,]
head(filtered_tpm_data)

# Extract genes that are common between the DESeq2 results and TPM data
matching_genes <- volcano_data[rownames(volcano_data) %in% rownames(filtered_tpm_data),]
head(matching_genes)

### Assign gene groups based on thresholds
# Create a new data frame with an initial group assignment for each gene
group_assignment <- as.numeric(matrix(1:nrow(matching_genes), nrow = nrow(matching_genes), ncol = 1))
gene_data <- as.data.frame(cbind(group = group_assignment, matching_genes))

# Set threshold values
lfc_threshold <- 1.0  # Threshold for log2 Fold Change
p_adjust_threshold <- 0.05  # p.adjust significance threshold

# Assign group labels based on log2FoldChange and p.adjust values
for (i in 1:nrow(gene_data)) {
    if (gene_data[i, 2] < -log10(p_adjust_threshold)) {  # If p.adjust > 0.05
        gene_data[i, 1] <- "A"  # Assign to group A (not significant)
    } else {
        if (gene_data[i, 3] > lfc_threshold) {  # If log2FoldChange > 1.0
            gene_data[i, 1] <- "B"  # Assign to group B (upregulated)
        } else if (gene_data[i, 3] < -lfc_threshold) {  # If log2FoldChange < -1.0
            gene_data[i, 1] <- "C"  # Assign to group C (downregulated)
        } else {
            gene_data[i, 1] <- "A"  # Assign to group A (not significant)
        }
    }
}

# Check the number of genes in each group
nrow(gene_data[gene_data$group == "A", ])  # Non-significant genes
nrow(gene_data[gene_data$group == "B", ])  # Upregulated genes
nrow(gene_data[gene_data$group == "C", ])  # Downregulated genes

### Create Volcano Plot
# Load required libraries
library(ggplot2)
library(reshape2)
library(ggrepel)

# Define labels for x and y axes
x_label <- expression(paste("", {Log[2]}, " Fold Change", sep = ""))  # X-axis label
y_label <- expression(paste("-", {Log[10]}, " (p.adjust)", sep = ""))  # Y-axis label

# Generate the plot
volcano_plot <- ggplot(gene_data, aes(x = log2FoldChange, y = padj, color = group, fill = group)) + 
        geom_point(size = 3.1, shape = 21) +  # Plot points with specified size and shape
        ggtitle("Volcano Plot: Case vs Control") +  # Add plot title
        labs(x = x_label, y = y_label) +  # Add x and y axis labels
        scale_x_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10)) +  # Set x-axis range and ticks
        scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +  # Set y-axis range and ticks
        theme(  # Customize plot appearance
                plot.title = element_text(family = "Arial", size = 12),  # Title font size
                axis.text = element_text(family = "Arial", size = 12, colour = "black"),  # Axis labels
                axis.title = element_text(family = "Arial", size = 12),  # Axis titles
                panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Black border around plot area
                axis.ticks = element_line(colour = "black", size = 1),  # Style of axis ticks
                legend.position = "none"  # Remove legend
        ) +
        geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dashed", size = 0.5) +  # Horizontal line at p.adjust threshold
        geom_vline(xintercept = lfc_threshold, colour = "black", linetype = "dashed", size = 0.5) +  # Vertical line at positive fold change threshold
        geom_vline(xintercept = -lfc_threshold, colour = "black", linetype = "dashed", size = 0.5) +  # Vertical line at negative fold change threshold
        annotate("rect", xmin = lfc_threshold, xmax = Inf, ymin = -log10(0.05), ymax = Inf, alpha = 0.15) +  # Highlight upregulated region
        annotate("rect", xmin = -lfc_threshold, xmax = -Inf, ymin = -log10(0.05), ymax = Inf, alpha = 0.15) +  # Highlight downregulated region
        scale_colour_manual(values = c("A" = "gray29", "B" = "darkred", "C" = "blue")) +  # Assign colors to groups
        scale_fill_manual(values = c("A" = "gray54", "B" = "red", "C" = "dodgerblue2"))

# Display the plot
volcano_plot

### Save the Volcano Plot
# Set the directory where the plot will be saved
output_directory <- "path/to/your/working/directory"  # Replace with your directory path
setwd(output_directory)

# Save the plot as a PNG file
ggsave(filename = "VolcanoPlot_Case_vs_Control.png", 
       plot = volcano_plot, 
       device = "png", 
       scale = 1, 
       width = 4, height = 4.5, 
       units = "in", 
       dpi = 900)
