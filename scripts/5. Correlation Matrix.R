### CorrelationMatrix

# 今までのデータ削除
rm(list=ls())

##################### リードデータの読み込み
path <- "path/to/your/working/directory"
setwd(path)


##
FileName01 <- "DataMatrix_ImportSampleList_RE-_vs_RE+_20240821AT.txt"
List0 <- read.table(FileName01, header=T, sep="\t", stringsAsFactors=F)
List1 <- List0[,c(1:10)] # abundanceのみを抽出する
colnames(List1)
head(List1,3)

FileName02 <- "Import_Sample_List_RE-_vs_RE+_20240821AT.txt"
List02 <- read.table(FileName02, header=T, sep="\t", stringsAsFactors=F)
colnames(List1) <- List02$sample
head(List1,3)
dim(List1)

TPM <- List1

# TPMが全てのサンプルいおいて1以上の発現を満たす遺伝子を抽出
library(genefilter)
f1 <- kOverA(10, A=1)        #「11組のサンプル以上でTPM1以上を持つ遺伝子を抽出」という条件（filter）をf1に格納
ffun <- filterfun(f1)        # フィルタリング用の関数(filtering function)を作成しffunに格納
obj <- genefilter(TPM, ffun) # 条件を満たすかどうかを判定した結果をobjに格納
obj                          # 確認
TPM.over1 <- TPM[obj,]       # objがTRUEとなる要素のみ抽出した結果をdataに格納
dim(TPM.over1)               # 確認

data1 <- TPM.over1
data2 <- log10(data1 + 1)
dim(data2)


# Calculate correlation coefficients
data3 <- cor(data2, method="spearman")
round(data3, digits=2)


library(corrplot)
packageVersion("corrplot")

col1 <- colorRampPalette(c("blue", "cyan", "white", "yellow", "red"))



corrplot(data3,
         method="color", # shade
         shade.col=NA,
         tl.col="black",
         
         number.cex = 1.5,
         addCoef.col="black",
         
         tl.cex = 1.2,    
         cl.cex = 1.2,  
         order = "hclust",
         col = col1(200),
         tl.srt=45,
         #mar = c(0, 0, 0, 0),
         col.lim = c(-1,1) #カラースケールの範囲
        )
