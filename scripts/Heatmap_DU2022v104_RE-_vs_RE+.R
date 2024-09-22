### Human RNA-seq

# Clear all previous data
rm(list=ls())

##################### Load read data
## Load the database
#DU2022
path <- "/Users/shara/Library/CloudStorage/DESeq2/DU2022_QC4_v104/gene_RE-_vs_RE+_20240829AT"
setwd(path)

Filename01　<- "DataMatrix_ImportSampleList_RE-_vs_RE+_20240821AT.txt"
list <- read.table(Filename01, header=T, sep="\t", stringsAsFactors=F)
dim(list)
head(list)
abundance <- list[,1:10]
head(abundance,3)

Filename02　<- "Import_Sample_List_RE-_vs_RE+_20240821AT.txt"
SAMPLE_list <- read.table(Filename02, header=T, sep="\t", stringsAsFactors=F)
colnames(abundance) <- SAMPLE_list$sample
head(abundance,3)


# Extract genes that meet the condition of TPM being 1 or higher in all samples
library(genefilter)
packageVersion("genefilter")
f1 <- kOverA(10, A=1)       
ffun <- filterfun(f1)        
obj <- genefilter(abundance, ffun)
obj           
abundance.over1 <- abundance[obj,]  
dim(abundance.over1)               


# Extract expression levels of the relevant genes
data1 <- abundance.over1
dim(data1)
head(data1)
data2 <- log10(data1 + 1)
dim(data2)


### Create Heatmap
x <- as.matrix(data2)

dev.off() # Close previous windows

d1 <- dist(x, method="manhattan")
d2 <- dist(t(x), method="manhattan")
# Select similarity or distance: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"

c1 <- hclust(d1, method="ward.D2")
c2 <- hclust(d2, method="ward.D2")
# Algorithm for clustering: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC), "centroid" (= UPGMC)

#install.packages("gplots") #Install
library(gplots)            #Load
library(RColorBrewer)
packageVersion("RColorBrewer")
packageVersion("gplots")


#Red and Blue, ward.D2
heatmap.2(x,
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),
          scale="row",
          key=T,
          keysize=1.0,
          density.info="none",
          trace="none",
          cexCol=1,
          cexRow=0.3,
          margins=c(3,10),
          Rowv=as.dendrogram(c1),
          Colv=as.dendrogram(c2)
)
