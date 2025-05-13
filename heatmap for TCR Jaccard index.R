rm(list = ls())
gc()
library(dplyr)
library(tidyr)
library(pheatmap)
LNCD8T <- readRDS("~/R/ICBstudy/result/data/Tcell/LNCD8T.rds")
tumorCD8Tcell <- readRDS("~/R/ICBstudy/result/data/Tcell/tumorCD8Tcell.rds")
LN<-LNCD8T
tumor<-tumorCD8Tcell
table(LN$TCR)
table(tumor$TCR)
LN<-LN[,LN$TCR%in%c("A","AAB","AB","B")]
tumor<-tumor[,tumor$TCR%in%c("A","AAB","AB","B")]
LN<-LN[,LN$Celltype!="MAIT"]
tumor<-tumor[,tumor$Celltype!="MAIT"]
LN_metadata <- LN@meta.data
Tumor_metadata <-tumor@meta.data



LN_tcr_counts <- LN_metadata %>%
  group_by(Celltype, CTnt) %>%
  summarise(count = n(), .groups = "drop")

Tumor_tcr_counts <- Tumor_metadata %>%
  group_by(Celltype, CTnt) %>%
  summarise(count = n(), .groups = "drop")


ln_celltypes <- unique(LN_metadata$Celltype)
tumor_celltypes <- unique(Tumor_metadata$Celltype)

similarity_matrix <- matrix(0, nrow = length(ln_celltypes), ncol = length(tumor_celltypes))
rownames(similarity_matrix) <- ln_celltypes
colnames(similarity_matrix) <- tumor_celltypes


for (ln_cell in ln_celltypes) {
  for (tumor_cell in tumor_celltypes) {
    ln_tcr <- LN_tcr_counts %>% filter(Celltype == ln_cell)
    tumor_tcr <- Tumor_tcr_counts %>% filter(Celltype == tumor_cell)
    common_tcr <- intersect(ln_tcr$CTnt, tumor_tcr$CTnt)
    common_count <- length(common_tcr)
    total_union <- length(union(ln_tcr$CTnt, tumor_tcr$CTnt))
    similarity <- ifelse(total_union == 0, 0, common_count / total_union)
    similarity_matrix[ln_cell, tumor_cell] <- similarity
  }
}
similarity_matrix1<-similarity_matrix[,c(2,1,3,4,5)]








p<-pheatmap(similarity_matrix1,
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(rev(c("#F39865", "white", "#7Fabd1")))(100),
         main = "TCR Jaccard Index for CD8+ Tcell Between TDLN and Tumor",
         fontsize_row = 15,
         fontsize_col = 15,
         cellwidth = 75,
         cellheight = 75,
         fontfamily = "Arial")

library(grid)
tiff("heatmap.tiff", width = 20, height = 18, units = "cm", 
     res = 1800, compression = "lzw", bg = "white")
grid.draw(p$gtable)  
dev.off()  




