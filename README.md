# Transcriptional profiling of patients with Fuchs endothelial corneal dystrophy with and without trinucleotide repeat expansion in TCF4


## Overview

This repository contains the analysis pipeline for transcriptional profiling of corneal endothelial cells (CECs) from patients with Fuchs Endothelial Corneal Dystrophy (FECD). The project compares transcriptional profiles between patients with and without trinucleotide repeat expansions in the TCF4 gene, using RNA sequencing (RNA-Seq) data.

This analysis uses various R packages, such as **DESeq2** and **tximport**, for differential gene expression analysis and visualization of the results, including volcano plots, PCA plots, heatmaps, and correlation matrices.

---

## Contents

| File | Description |
|------|-------------|
| `SRRXXXXXXXX.genes.results` | Gene expression results from RNA-Seq experiments, generated by RSEM, which estimates gene and isoform expression levels. |
| `Sample_List_YYYYMMDD.txt` | Metadata file containing information about the samples used in the RNA-Seq experiment. |
| `DataMatrix_ImportSampleList_Case_vs_Control.txt` | Gene expression data matrix, containing expression levels (e.g., TPM) for each gene in each sample. |
| `DESeq2_Result_of_WaldTest_Case_vs_Control.txt` | Results from differential expression analysis performed using the DESeq2 package, comparing case samples to control samples. |
| 1. `DESeq2` | Script for performing differential expression analysis using DESeq2 and generating results with log2 fold changes and p-values. |
| 2. `VolcanoPlot.R` | Script for generating a volcano plot from differential expression results. |
| 3. `PCA.R` | Script for generating a PCA plot of the gene expression data. |
| 4. `Heatmap.R` | Script for creating a heatmap of differentially expressed genes. |
| 5. `CorrelationMatrix.R` | Script for generating a correlation matrix of the gene expression data. |
| 6. `GO_enrichment_analysis.R` | Script for performing Gene Ontology (GO) enrichment analysis using the differentially expressed genes. |

---

### Example of Table Format of Sample_List_YYYYMMDD.txt
"SRRXXXXXXXX.genes.results" files were obtained by RSEM

|  sample  | group |  RIN  | path |
|----------|-------|-------|------|
| Case1  | case  | 6.3   | /path/to/your/working/directory/SRRXXXXXXXX.genes.results |
| Case2  | case  | 8.0   | /path/to/your/working/directory/SRRXXXXXXXX.genes.results |
| Case3  | case  | 7.3   | /path/to/your/working/directory/SRRXXXXXXXX.genes.results |
| Control1 | ctr | 7.4   | /path/to/your/working/directory/SRRXXXXXXXX.genes.results |
| Control2 | ctr | 8.0   | /path/to/your/working/directory/SRRXXXXXXXX.genes.results |
| Control3 | ctr | 7.2   | /path/to/your/working/directory/SRRXXXXXXXX.genes.results |

---

## Requirements

### Software
- **R** version 4.3.2 or later
- **R Packages**:
  - `DESeq2`
  - `tximport`
  - `ggplot2` (for plotting)
  - `genefilter` (for filtering genes)
  - `corrplot` (for correlation matrix plots)
  - `gplots` and `RColorBrewer` (for heatmaps)

---

### Data
- RNA-Seq gene expression data in `SRRXXXXXXXX.genes.results` format (generated by RSEM or similar tools).
- A sample list file (`Sample_List_YYYYMMDD.txt`) containing metadata for all samples.

---

### Acquisition of RNA-Seq Data
The raw RNA-Seq FASTQ files used for this analysis were obtained from the DNA Data Bank of Japan (DDBJ) Sequence Read Archive (SRA). The data can be downloaded under the following accession IDs:
- Non-FECD (Control): DRP006678
- FECD Patients: DRA015078

---

## Citation

If you use this code or data, please cite our paper:

```
[Citation information will be updated soon]
```

---

## Contact

For questions or issues, please contact:

**Tatsuya Nakagawa**  
Email: tatsuya.n.nakagawa@gmail.com
