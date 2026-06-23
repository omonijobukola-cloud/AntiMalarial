#install deseq2

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")




#Load deseq2 library

library(DESeq2)
library(ggplot2)


#load the counts table

Counts <- read.delim("GSE198699_gene_count.csv", header = TRUE, row.names = 1, sep = ",")



#Filter the unnecessary data by using those greater than 50

Counts <- Counts[which(rowSums(Counts) > 50),]


#Call counts

Counts


#set condition

condition <- factor(c("A","A","A", "B","B","B", "C","C","C", "D","D"))


#make a data frame for conditions and samples

coldata <- data.frame(row.names = colnames(Counts), condition)

coldata



#start calling desyq functions

dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

dds <- DESeq(dds)

vsdata <- vst(dds, blind=FALSE)

plotPCA(vsdata, intgroup = "condition")


#check my data dispersion

plotDispEsts(dds)


#create a table of comparison for for the degs

res <- results(dds, contrast = c("condition", "A", "B"))

res


#filter insignificant genes

sigs <- na.omit(res)

sigs <- sigs[sigs$padj < 0.05,]

sigs

#save sigs as csv file
write.csv(sigs,file = "deseqresults.csv")


#Convert ensemble id to gene symbol

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("AnnotationDbi")

library(org.Mm.eg.db)
library(AnnotationDbi)

sigs.df <- as.data.frame(sigs)

sigs.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL" , column = "SYMBOL")

sigs.df


#PROCESS OF MAKING HEATMAP
#BiocManager::install("ComplexHeatmap")

sigs.df <- sigs.df[(sigs.df$baseMean > 450) & (abs(sigs.df$log2FoldChange) >3.5),]

library(ComplexHeatmap)

mat <- counts(dds,normalized = T)[rownames(sigs.df),]

mat.z <- t(apply(mat, 1, scale))
coldata

colnames(mat.z) <- rownames(coldata)
mat.z

#install.packages("magick")
#library(magick)

h <- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), name = "Z-score", row_labels = sigs.df[rownames(mat.z),]$symbol)

png('simple_heatmap.png', res = 250, width = 2000, height = 4000)
print(h)
dev.off()


#Making a volcano plot
BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)
EnhancedVolcano(sigs.df, x = "log2FoldChange", y = "padj", lab = sigs.df$symbol)
