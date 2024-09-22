### Created on: 20240829
### Volcano Plot

rm(list=ls()) 
#dev.off() # Close the window if there is any issue


##################### Read data
## Read Database
path <- "/Users/shara/Library/CloudStorage/OneDrive-同志社大学/遺伝子解析班データファイル/RNA-Seqプロジェクト/DESeq2/DU2022_QC4_v104/gene_RE-_vs_Control_20240829AT"
setwd(path) 


# Result of DESeq2
FileName <- "DESeq2_Result_of_WaldTest_RE-_vs_Control_20240829AT.txt"
data1 <- read.table(FileName, header=T, sep="\t", stringsAsFactors=F)
data2 <- cbind(padj = data1$padj, log2FoldChange = data1$log2FoldChange)
rownames(data2) <- rownames(data1)
dim(data2)
head(data2)

data3 <- data2
data3[,1] <- -log10(data2[,1])
head(data3,3)


# Expression level data
path.tpm <- "/Users/shara/Library/CloudStorage/OneDrive-同志社大学/遺伝子解析班データファイル/RNA-Seqプロジェクト/DESeq2/DU2022_QC4_v104/gene_FECD_vs_Control_20240912AT"
setwd(path.tpm)

FileName.tpm <- "DataMatrix_ImportSampleList_FECD_vs_Control_20240912AT.txt"  #ファイル名を入力
data1.tpm <- read.table(FileName.tpm, header=T, sep="\t", stringsAsFactors=F)
head(data1.tpm,3)

data2.tpm <- data1.tpm[,1:17]
dim(data2.tpm)
head(data2.tpm,3)

data3.tpm <- data2.tpm[rowMeans(data2.tpm) >= 1,]
dim(data3.tpm)
head(data3.tpm,3)

# Extract only matching gene names
data.list <- data3[rownames(data3) %in% rownames(data3.tpm),]
dim(data.list)


############################################################
#### Apply conditions based on LFC and P.adjust
nrow(data.list)

group <- as.numeric(matrix(1:nrow(data.list), nrow = nrow(data.list), ncol=1))
data4 <- as.data.frame(cbind(group,data.list))
dim(data4)

lfc.threshold <- 1.0 # Log2_Fold_Change threshold

x <- 1
while (x <= nrow(data.list)) {   # while ( condition )

        if (data4[x, 2] < -log10(0.05)) {  # If p.adjust is greater than 0.05
                data4[x, 1] <- "A"          # Assign A
        } else {
                if (data4[x, 3] > lfc.threshold) { 　# If Log2_Fold_Change is greater than 1.0
                        data4[x, 1] <- "B"  # Assign B
                } else {
                        if (data4[x, 3] < -lfc.threshold) { 　# If Log2_Fold_Change is less than -1.0
                                data4[x, 1] <- "C"    # Assign C
                        } else {
                                data4[x, 1] <- "A"   # Assign 1 if Log2_Fold_Change is between -1.0 and 1.0
                        }
                }
        }
  
        x <- x + 1 
}

# Check differentially expressed genes
nrow(data4[data4$group == "A", ])
nrow(data4[data4$group == "B", ])
nrow(data4[data4$group == "C", ])

######## Drawing the Volcano Plot #################################
library(ggplot2) # library
library(reshape2)
library(ggrepel)
packageVersion("ggplot2")
packageVersion("reshape2")
packageVersion("ggrepel")


x.lab = expression(paste("",{Log[2]}," ","Fold Change", sep="")) #x軸ラベル
y.lab = expression(paste("-",{Log[10]}," ","(p.adjust)", sep="")) #y軸ラベル


max(data4$padj)
min(data4$padj)
max(data4$log2FoldChange)
min(data4$log2FoldChange)


plot <- ggplot(data4, aes(x = log2FoldChange, y = padj, color = group, fill = group)) + # Define data frame to be used for plotting; define data for x and y axes; crate a scatterplot object.
        #geom_point() +
        geom_point(size = 3.1, shape = 21) + # Define data point style.
        
        ggtitle("DU2022v104 RE- vs Control") + # Define title and subtitle.
        
        labs(x = x.lab, y = y.lab) + # Define labels for x and y axes.
        
        scale_x_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10)) + # Define x limits, add ticks.
        scale_y_continuous(limits = c(0,15), breaks = seq(0, 15, by = 5)) + # Define y limits, add ticks.
        
        theme(
                plot.title = element_text(family = "Arial", size = 12, hjust = 0), # Title size and font.
                plot.subtitle = element_text(family = "Arial", size = 12), # Subtitle size and font.
                axis.text = element_text(family = "Arial", size = 12, colour = "black"), # Size and font of x and y values.
                axis.title = element_text(family = "Arial", size = 12), # Size and font of x and y axes.
                panel.border = element_rect(colour = "black", fill = NA, size = 1), # Black border around the plot area.
                axis.ticks = element_line(colour = "black", linewidth = 1), # Style of x and y ticks.
                legend.position = "none"
        ) + # Remove legend.
        
        geom_hline(yintercept = 1.30103, colour = "black", linetype = "dashed", size = 0.5) + # Horizontal significance cut-off line.
        geom_vline(xintercept = lfc.threshold, colour = "black", linetype = "dashed", size = 0.5) + # Vertical significance cut-off line (+).
        geom_vline (xintercept = -lfc.threshold, colour = "black", linetype = "dashed", size = 0.5) + #Vertical significance cut-off line (-)
        
        annotate("rect", xmin = lfc.threshold, xmax = Inf, ymin = 1.30103, ymax = Inf, alpha = .15) + #alpha (default: 1=opaque) transparency of the rectangle's fill
        annotate("rect", xmin = -lfc.threshold, xmax = -Inf, ymin = 1.30103, ymax = Inf, alpha = .15) +
        #annotate("text",family = "Arial", x = 3.65, y = 1.45, label = "P-value = 0.05", size = 4, fontface = "plain") + # Label to horizontal cut-off line.
        #annotate("text",family = "Arial", x = 0.39, y = 6.4, label = expression(paste("",{Log[2]}," ","Fold Change = 0.5", sep="")),
                 #size = 4, fontface = "plain", srt = 90) +# Label to vertical cut-off line.
        #annotate("text",family = "Arial", x = -0.6, y = 6.4, label = expression(paste("",{Log[2]}," ","Fold Change = -0.5", sep="")),
                 #size = 4, fontface = "plain", srt = 90) +# Label to vertical cut-off line.

        scale_colour_manual(values = c("A" = "gray29","B" = "darkred","C" = "blue")) +
        scale_fill_manual(values = c("A" = "gray54","B" = "red","C" = "dodgerblue2"))

plot

setwd("/Users/shara/Library/CloudStorage/OneDrive-同志社大学/遺伝子解析班データファイル/RNA-Seqプロジェクト/Noexpansion vs Expansion/論文用 RE- vs+/v104/DU2022/ボルケーノプロット/新")

ggsave(filename = "DU2022 RE- vs Control.png", 
       plot = plot, 
       device = "png", 
       scale = 1, 
       width = 4, height = 4.5, 
       units = c("in"),
       dpi = 900)

