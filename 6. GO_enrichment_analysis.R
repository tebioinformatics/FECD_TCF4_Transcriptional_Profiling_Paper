# GO enrichment analysis

# Clear all previous data and close any open graphical devices
rm(list=ls()) # Remove all objects from the environment
dev.off()     # Close all graphical devices

# Load necessary libraries
library(biomaRt)         # For accessing the Ensembl biomart database
library(org.Hs.eg.db)    # Human gene annotation package
library(DOSE)            # For semantic similarity calculation among GO terms
library(ggplot2)         # For visualization
library(clusterProfiler) # For functional enrichment analysis
library(enrichplot)      # For enrichment plot visualization
library(forcats)         # For factor reordering
library(dplyr)           # Data manipulation
library(scales)          # For scaling functions
library(tidyverse)       # Collection of data science tools

# Set working directory and load DESeq2 results
working_directory <- "path/to/your/working/directory"
setwd(working_directory)

# Load DESeq2 results of Wald Test (differential expression analysis)
deseq2_filename <- "DESeq2_Result_of_WaldTest_Case_vs_Control.txt"
deseq2_data <- read.table(deseq2_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ensembl_gene_id <- rownames(deseq2_data)
deseq2_data <- data.frame(ensembl_gene_id, deseq2_data)

# Load expression level data (TPM) and filter genes with average expression >= 1
tpm_filepath <- "/Users/"
setwd(tpm_filepath)
tpm_filename <- "DataMatrix_YYYYMMDD.txt"
tpm_data <- read.table(tpm_filename, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
filtered_tpm_data <- tpm_data[rowMeans(tpm_data[,1:17]) >= 1,]

# Filter significant genes based on adjusted p-value and log2FoldChange threshold
significant_genes <- deseq2_data[rownames(deseq2_data) %in% rownames(filtered_tpm_data),]
significant_genes <- significant_genes[significant_genes$padj < 0.05,]
logfc_threshold <- 1.0
upregulated_genes <- significant_genes[significant_genes$log2FoldChange > logfc_threshold,]
downregulated_genes <- significant_genes[significant_genes$log2FoldChange < -logfc_threshold,]

# Extract Ensembl gene IDs for further analysis
upregulated_gene_ids <- upregulated_genes$ensembl_gene_id
downregulated_gene_ids <- downregulated_genes$ensembl_gene_id
all_regulated_gene_ids <- c(upregulated_gene_ids, downregulated_gene_ids)

# Use biomaRt to convert Ensembl gene IDs to Entrez IDs and other annotations
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
convert_genes <- function(gene_ids) {
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "entrezgene_description", 
                       "hgnc_symbol", "chromosome_name", "start_position", "end_position", 
                       "band", "gene_biotype"),
        filters = "ensembl_gene_id", values = gene_ids, mart = mart)
}

upregulated_converted <- convert_genes(upregulated_gene_ids)
downregulated_converted <- convert_genes(downregulated_gene_ids)
all_converted <- convert_genes(all_regulated_gene_ids)

# Filter for protein-coding genes and remove missing values
protein_coding_filter <- function(df) df[df$gene_biotype == "protein_coding",]
upregulated_protein_coding <- protein_coding_filter(upregulated_converted)
downregulated_protein_coding <- protein_coding_filter(downregulated_converted)
all_protein_coding <- protein_coding_filter(all_converted)

# GO enrichment analysis using clusterProfiler
go_enrich <- function(gene_ids) {
  enrichGO(gene_ids, OrgDb = org.Hs.eg.db, ont = "all", readable = TRUE, 
           pAdjustMethod = "BH", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
}

go_results_up <- go_enrich(upregulated_protein_coding$entrezgene_id)
go_results_down <- go_enrich(na.omit(downregulated_protein_coding$entrezgene_id))
go_results_all <- go_enrich(na.omit(all_protein_coding$entrezgene_id))

# Process and filter GO analysis results
go_results_df <- as.data.frame(go_results_all)
significant_go_results <- go_results_df[go_results_df$p.adjust < 0.05,]

# Extract specific terms for Biological Process (BP), Cellular Component (CC), and Molecular Function (MF)
extract_top_terms <- function(go_df, ontology, top_n) {
  ontology_df <- go_df[go_df$ONTOLOGY == ontology, ]
  sorted_ontology_df <- ontology_df[order(ontology_df$p.adjust, ontology_df$pvalue),]
  head(na.omit(sorted_ontology_df), top_n)
}

top_bp <- extract_top_terms(significant_go_results, "BP", 8)
top_cc <- extract_top_terms(significant_go_results, "CC", 8)
top_mf <- extract_top_terms(significant_go_results, "MF", 8)

# Combine top BP, CC, MF terms into one data frame
combine_terms <- function(category_name, terms_df) {
  category_col <- matrix(rep(category_name, nrow(terms_df)), ncol = 1)
  colnames(category_col) <- "Category"
  cbind(category_col, terms_df)
}

combined_go_results <- rbind(combine_terms("Biological Process", top_bp),
                             combine_terms("Cellular Component", top_cc),
                             combine_terms("Molecular Function", top_mf))

# Save the results to a file
output_directory <- "/path/to/your/output/directory"
setwd(output_directory)
write.table(combined_go_results, file = "GO_analysis_results.txt", sep = "\t", quote = FALSE)

# Visualization of GO enrichment results using ggplot2
plot <- ggplot(combined_go_results, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  facet_grid(Category ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(color = "black", size = 1.5),
        panel.grid = element_line(size = 0.1, linetype = "dotted"),
        strip.text = element_text(size = 13))

# Save the plot
ggsave(filename = "GO_enrichment_plot.png", plot = plot, device = "png", 
       width = 10, height = 8, dpi = 900)
