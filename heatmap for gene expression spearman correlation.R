
library(Seurat)
library(dplyr)
library(tidyr)
library(pheatmap)
LNCD8T <- readRDS("~/R/ICBstudy/result/data/Tcell/LNCD8T.rds")
tumorCD8Tcell <- readRDS("~/R/ICBstudy/result/data/Tcell/tumorCD8Tcell.rds")

tumor<-tumorCD8Tcell
LN<-LNCD8T
tumor<-tumor[,tumor$Celltype!="MAIT"]
LN<-LN[,LN$Celltype!="MAIT"]




LN_avg_expr <- AverageExpression(LN, group.by = "Celltype", assays = "RNA")$RNA
LN_avg_expr <- as.data.frame(LN_avg_expr)  


tumor_avg_expr <- AverageExpression(tumor, group.by = "Celltype", assays = "RNA")$RNA
tumor_avg_expr <- as.data.frame(tumor_avg_expr)  


common_genes <- intersect(rownames(LN_avg_expr), rownames(tumor_avg_expr))
LN_avg_expr <- LN_avg_expr[common_genes, ]
tumor_avg_expr <- tumor_avg_expr[common_genes, ]

LN_celltypes <- colnames(LN_avg_expr)
tumor_celltypes <- colnames(tumor_avg_expr)
correlation_matrix <- matrix(0, nrow = length(LN_celltypes), ncol = length(tumor_celltypes),
                             dimnames = list(LN_celltypes, tumor_celltypes))


for (LN_cell in LN_celltypes) {
  for (tumor_cell in tumor_celltypes) {
    correlation_matrix[LN_cell, tumor_cell] <- cor(LN_avg_expr[[LN_cell]], 
                                                   tumor_avg_expr[[tumor_cell]], 
                                                   method = "spearman")
  }
}

correlation_matrix<-t(correlation_matrix)


correlation_matrix1<-correlation_matrix[c(5,1,2,4,3),c(4,1,2,3)]

p<-pheatmap(correlation_matrix1,
            cluster_rows = F,
            cluster_cols = F,
            color = colorRampPalette(rev(c("#F39865", "white", "#7Fabd1")))(100),
            main = "Gene Expression Spearman Correlation",
            fontsize_row = 10,
            fontsize_col = 10,
            cellwidth = 30,
            cellheight = 30,
            fontfamily = "Arial")





library(grid)
tiff("genesimilarity.tiff", width = 20, height = 12, units = "cm", 
     res = 3600, compression = "lzw", bg = "white")
grid.draw(p$gtable) 
dev.off()  

