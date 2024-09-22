### Creation Date: 2024-08-29


rm(list=ls()) # Clear all previous data
############### Installing the Heatplus package (only for the first time) 
#install.packages("maptools")
#install.packages("scatterplot3d")

############### Font settings
par(family="sans")  # Set graph font to Times New Roman; serif is Times New Roman, sans is Arial
par(pty="s")  # Make the plot square-shaped

############### Loading read data
path <- "/Users/shara/Library/CloudStorage/DESeq2/DU2022_QC4_v104/gene_RE-_vs_RE+_20240829AT"
setwd(path)  # 作業ディレクトリを変更する

# TPM
FileName2 <- "DataMatrix_ImportSampleList_RE-_vs_RE+_20240821AT.txt"
Data_matrix0 <- read.table(FileName2, sep = "\t", header = T)
TPM <- Data_matrix0[,c(1:10)]
head(TPM)
dim(TPM)
colnames(TPM)

# Extract genes that meet the condition of TPM being 1 or higher in all samples
library(genefilter)
packageVersion("genefilter")
f1 <- kOverA(10, A=1)  # Set the condition (filter) to extract genes that have a TPM of 1 or more in over 11 samples, store it in f1
ffun <- filterfun(f1)  # Create a filtering function and store it in ffun
obj <- genefilter(TPM, ffun)  # Judge whether the condition is met, store the result in obj
obj  # Check
TPM.over1 <- TPM[obj,]  # Store the elements that are TRUE in obj in the data
dim(TPM.over1)  # Check


### Extract data from the data_matrix
TPM.DEGs <- TPM.over1
dim(TPM.DEGs)
head(TPM.DEGs)


FileName3 <- "Import_Sample_List_RE-_vs_RE+_20240821AT.txt"
sample_list<- read.table(FileName3, sep = "\t", header = T)

colnames(TPM.DEGs) <- sample_list$sample
head(TPM.DEGs)

############### PCA (2 groups) using the prcomp function
TPM.DEGs.log10 <- log10(TPM.DEGs + 1) # Log transformation

# PCA calculation
pcaobj <- prcomp(t(TPM.DEGs.log10),  # Transpose to make cohorts vertical and genes horizontal
                 scale = FALSE)  # Set scale to FALSE as the units are aligned


# Proportion of Variance (axis contribution)
summary(pcaobj)$importance

# PCA plot (two-dimensional) (PC1, 2)
PC1 <- pcaobj$x[, 1]
PC2 <- pcaobj$x[, 2]


x.1 <- data.frame(PC1 = c(PC1), PC2 = c(PC2))

Nex <- matrix(1:4, nrow=4, ncol=1)
Nex[,1] <- "RE-"
Ex <- matrix(1:6, nrow=6, ncol=1)
Ex[,1] <- "RE+"
Group <- rbind(Nex, Ex)

Sample <- rownames(x.1)
x.2 <- cbind(Sample, x.1, Group)
head(x.2)

# Axis labels for contribution rate
kiyo1 <- signif(summary(pcaobj)$importance[2,1]*100, 3)
kiyo2 <- signif(summary(pcaobj)$importance[2,2]*100, 3)

PC1.lab <- paste("PC1  ", "(", kiyo1, "%)", sep = "")
PC2.lab <- paste("PC2  ", "(", kiyo2, "%)", sep = "")

packageVersion("ggplot2")

# Plot the graph
library(ggplot2)
ggplot(x.2 , aes(x = PC1, y = PC2, group = Group, fill = Group))+
  theme_linedraw()+
  
  labs(x=PC1.lab)+
  labs(y=PC2.lab)+
  
  geom_text(aes(label = Sample), vjust = -1, size = 5)+

  geom_point(aes(shape = Group), size = 7, stroke = 0.8)+

  scale_shape_manual(values = c("RE-"=23,"RE+"=21))+
  scale_fill_manual(values = c("#67A9CF", "#EF8A62"))+
  
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1.0), # Black border around the plot area.
        panel.grid = element_line(size=1.0, colour = "black", linetype = "dashed"), # dotted, solid, blank, dashed
        aspect.ratio = 1.0, # アスペクト比)
        text = element_text(size = 17)
  )+
  
  xlim(-20,20)+
  ylim(-20,20)
