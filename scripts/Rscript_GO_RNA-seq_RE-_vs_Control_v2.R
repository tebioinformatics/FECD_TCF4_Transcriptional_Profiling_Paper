#### Created on: 20240130
### GO analysis
### RE- vs Control

rm(list=ls()) # Clear all previous data
dev.off() # Close windows if there are any issues

##################### Load the data
## Load the database
path <- "/Users/shara/Library/CloudStorage/OneDrive-同志社大学/遺伝子解析班データファイル/RNA-Seqプロジェクト/DESeq2/DU2022_QC4_v104/gene_RE-_vs_Control_20240829AT"

setwd(path)      # Change the working directory
FileName <- "DESeq2_Result_of_WaldTest_RE-_vs_Control_20240829AT.txt"  # Specify the file name

# Load the file data into data0
data1 <- read.table(FileName, header=T, sep="\t", stringsAsFactors=F)

head(data1)
dim(data1)

ensembl_gene_id <- rownames(data1) # Get rownames
data2 <- data.frame(ensembl_gene_id, data1)
head(data2)
dim(data2)

########## Extract differentially expressed genes ############
# Expression level data
path.tpm <- "/Users/shara/Library/CloudStorage/OneDrive-同志社大学/遺伝子解析班データファイル/RNA-Seqプロジェクト/DESeq2/DU2022_QC4_v104/gene_FECD_vs_Control_20240912AT"
setwd(path.tpm)

FileName.tpm <- "DataMatrix_ImportSampleList_FECD_vs_Control_20240912AT.txt"  
data1.tpm <- read.table(FileName.tpm, header=T, sep="\t", stringsAsFactors=F)
head(data1.tpm,3)

data2.tpm <- data1.tpm[,1:17]
dim(data2.tpm)
head(data2.tpm,3)

data3.tpm <- data2.tpm[rowMeans(data2.tpm) >= 1,] # Extract genes with average expression level of 1 or higher
dim(data3.tpm)
head(data3.tpm,3)


data_expression <- data2[rownames(data2) %in% rownames(data3.tpm),] # Extract genes with an average expression level of 1 or higher & remove missing values
head(data_expression,3)
dim(data_expression)


data_pvalue <- data_expression[data_expression$padj < 0.05, ] # Extract significant genes
dim(data_pvalue)
head(data_pvalue)

th.lfc <- 1.0 # Set threshold
upregulated   <- data_pvalue[data_pvalue$log2FoldChange > th.lfc, ]
downregulated <- data_pvalue[data_pvalue$log2FoldChange < -th.lfc, ]  
dim(upregulated)
head(upregulated)
dim(downregulated)


upregulated.genes   <- c(upregulated$ensembl_gene_id) # Extract genes (ENTREZ) for analysis
downregulated.genes <- c(downregulated$ensembl_gene_id)
regulated.genes <- c(upregulated.genes,downregulated.genes)
length(regulated.genes)

non_regulated   <- data2[!(data2$padj < 0.05 & abs(data2$log2FoldChange) > th.lfc), ]
dim(non_regulated)


# Convert and obtain gene IDs
#BiocManager::install("biomaRt")
#packageVersion("biomaRt")
library(biomaRt)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# upregulated
up.conv <- getBM(attributes = c("ensembl_gene_id","entrezgene_id","entrezgene_description","hgnc_symbol","chromosome_name","start_position","end_position","band","gene_biotype"),
                 filters    = "ensembl_gene_id",
                 values     = upregulated.genes, 
                 mart       = mart)
dim(up.conv)
head(up.conv)

# downregulated
down.conv <- getBM(attributes = c("ensembl_gene_id","entrezgene_id","entrezgene_description","hgnc_symbol","chromosome_name","start_position","end_position","band","gene_biotype"),
                   filters    = "ensembl_gene_id",
                   values     = downregulated.genes, 
                   mart       = mart)
dim(down.conv)
head(down.conv)

all.conv <- getBM(attributes = c("ensembl_gene_id","entrezgene_id","entrezgene_description","hgnc_symbol","chromosome_name","start_position","end_position","band","gene_biotype"),
                   filters    = "ensembl_gene_id",
                   values     = regulated.genes, 
                   mart       = mart)
dim(all.conv)
head(all.conv)


# Extract only protein_coding genes
non_protein.up <- up.conv[up.conv$gene_biotype != "protein_coding",]
dim(non_protein.up)
protein.up <- up.conv[up.conv$gene_biotype == "protein_coding",]
dim(protein.up)

DEGs1.up.exNA <- na.omit(protein.up$entrezgene_id)
length(DEGs1.up.exNA)


non_protein.down <- down.conv[down.conv$gene_biotype != "protein_coding",]
dim(non_protein.down)
protein.down <- down.conv[down.conv$gene_biotype == "protein_coding",]
dim(protein.down)
DEGs1.down.exNA <- na.omit(protein.down$entrezgene_id)
length(DEGs1.down.exNA)

non_protein.all <- all.conv[all.conv$gene_biotype != "protein_coding",]
dim(non_protein.all)
protein.all <- all.conv[all.conv$gene_biotype == "protein_coding",]
dim(protein.all)
DEGs1.all.exNA <- na.omit(protein.all$entrezgene_id)
length(DEGs1.all.exNA)



### install library
#BiocManager::install("ggplot2")
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
packageVersion("ggplot2")
packageVersion("org.Hs.eg.db")
packageVersion("DOSE")

#BiocManager::install("clusterProfiler") # install
library(clusterProfiler)
packageVersion("clusterProfiler")

#BiocManager::install("enrichplot") # for gseaplot2
library(enrichplot)
packageVersion("enrichplot")

### GO analysis (enrichGO)
go.up <- enrichGO(c(protein.up$entrezgene_id),
                  OrgDb = "org.Hs.eg.db",
                  ont="all",   # all, BP, MF, CC
                  readable = TRUE,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.5,
                  qvalueCutoff  = 0.5)

go.down <- enrichGO(DEGs1.down.exNA,
                    OrgDb = "org.Hs.eg.db",
                    ont="all",
                    readable = TRUE,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.5,
                    qvalueCutoff  = 0.5)

go.all <- enrichGO(DEGs1.all.exNA,
                    OrgDb = "org.Hs.eg.db",
                    ont="all",
                    readable = TRUE,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.5,
                    qvalueCutoff  = 0.5)




### Process the analysis results

box0   <- as.data.frame(go.all)
dim(box0)

box1   <- box0[box0$p.adjust < 0.05,]
dim(box1)


num.BP.1 <- 8   # Number of terms to display for BP
num.CC.1 <- 8
num.MF.1 <- 8


### Extract BP, CC, MF
BP <- box1[box1$ONTOLOGY == "BP", ]
BP01 <- BP[order(BP$p.adjust,BP$pvalue, decreasing = c(FALSE,FALSE)),]
BP02 <- na.omit(BP01[c(1:num.BP.1),]); #必要な列を抽出
dim(BP)
head(BP02)

CC <- box1[box1$ONTOLOGY == "CC", ]
CC01 <- CC[order(CC$p.adjust,CC$pvalue, decreasing = c(FALSE,FALSE)),]
CC02 <- na.omit(CC01[c(1:num.CC.1),]); #必要な列を抽出
dim(CC)

MF <- box1[box1$ONTOLOGY == "MF", ]
MF01 <- MF[order(MF$p.adjust,MF$pvalue, decreasing = c(FALSE,FALSE)),]
MF02 <- na.omit(MF01[c(1:num.MF.1),]); #必要な列を抽出
dim(MF)

# Create category boxes
categ.BP <- matrix(1:num.BP.1, nrow=num.BP.1, ncol=1) #1×nの行列を作製する
categ.BP[,1] <- "Biological Process"
colnames(categ.BP) <- c("Category")

categ.CC <- matrix(1:num.CC.1, nrow=num.CC.1, ncol=1) #1×nの行列を作製する
categ.CC[,1] <- "Cellular Component"
colnames(categ.CC) <- c("Category")

categ.MF <- matrix(1:num.MF.1, nrow=num.MF.1, ncol=1) #1×nの行列を作製する
categ.MF[,1] <- "Molecular Function"
colnames(categ.MF) <- c("Category")


# Combine the data
BP03   <- cbind(categ.BP, BP02)
CC03   <- cbind(categ.CC, CC02)
MF03   <- cbind(categ.MF, MF02)


GO.result <- as.data.frame(rbind(BP03, CC03, MF03))
head(GO.result)



#### Caluculate GeneRatio
#install.packages("tidyverse")
packageVersion("tidyverse")
library(tidyverse)
GeneRatio0 <- as.matrix(str_split_fixed(GO.result$GeneRatio, pattern = "/", n = 2))　# Split by "/"
num.2 = nrow(GeneRatio0)

GeneRatio.matrix <- matrix(1:num.2, nrow=num.2, ncol=1) #1×nの行列を作製する

i =1

for(i in 1:num.2) {
  GeneRatio.matrix[i,1] <- (as.numeric(GeneRatio0[i,1])/as.numeric(GeneRatio0[i,2]))
}

GO.result01 <- as.data.frame(cbind(GO.result, GeneRatio.matrix))


#### Convert text to uppercase
box2 <- as.matrix(str_split_fixed(GO.result01$Description, pattern = "", n = 3))

#install.packages("stringr")
packageVersion("stringr")
library(stringr)
num.3 = nrow(box2)

for(i in 1:num.3) {
  if (box2[i,2] == "R") { # If "R" is part of "RNA", do nothing
    
  } else {
    box2[i,1] <- toupper(box2[i,1]) # Convert to uppercase
  }
} 

# Combine text
for(i in 1:num.3) {
  box2[i,1] <- paste(box2[i,1], box2[i,2], sep = "")
  box2[i,1] <- paste(box2[i,1], box2[i,3], sep = "")
}


GO.result01$Description <- box2[,1]
head(GO.result01)

# Log-transform p.adjust
GO.result02 <- GO.result01
GO.result02$p.adjust <- -log10(GO.result01$p.adjust)


# Extract specific GO terms
GO.result03 <- GO.result02


setwd("/Users/shara/Library/CloudStorage/OneDrive-同志社大学/遺伝子解析班データファイル/RNA-Seqプロジェクト/Noexpansion vs Expansion/論文用 RE- vs+/v104/DU2022/GO解析/RIN/RE- vs Control")
write.table(GO.result03, file="GO_analysis_downregulated_genes_RE-_vs_Control.txt", sep="\t", quote=F)

# load the library
library(forcats)
library(dplyr)
library(scales)
packageVersion("forcats")
packageVersion("dplyr")
packageVersion("scales")



### Visualization
plot <- ggplot(GO.result03, aes(x=p.adjust , y= reorder(Description,p.adjust),fill=p.adjust))+  # Reverse order by placing "-" before GeneRatio
  #geom_point()+
  geom_bar(stat = "identity", width = 0.6, color = "black", fill="#67A9CF",size = 0.6)+ # Specify bar width
  
  facet_grid(Category~., scale="free")+
  
  theme_bw()+ # white background and gray grid lines，テーマを先に指定してからフォントと色を変更する
  
  labs(x="")+
  labs(y="")+
  
  theme(  axis.text.x = element_text(family = "Arial", color = "black", size=15),  # Axis label font size
          axis.text.y = element_text(family = "Arial", color = "black", size=10), 
          axis.title = element_text(family = "Arial", color = "black", size=10), # Axis text font size
          legend.title = element_text(family = "Arial", color = "black", size=10),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5), # Black border around the plot area.
          panel.grid = element_line(size=0.1, colour = "black", linetype = "dotted"), # dotted, solid, blank, dashed
          strip.text = element_text(family = "Arial", color = "black", size = 13),
          aspect.ratio = 1.8 # アスペクト比
  )
plot 

ggsave(filename = "DU2022_RE-_vs_Control_RIN_ver_Down.png", 
       plot = plot, 
       device = "png", 
       scale = 1, 
       width = 10, height = 8, 
       units = c("in"),
       dpi = 900)

