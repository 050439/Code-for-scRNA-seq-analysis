rm(list = ls())
gc()
tumorCD8Tcell <- readRDS("~/R/ICBstudy/result/data/Tcell/tumorCD8Tcell.rds")
table(tumorCD8Tcell$groupdetail)
tumorCD8Tcell$groupdetail[tumorCD8Tcell$groupdetail=="dMMR-naive"]<-"dMMR_ICI(-)"
tumorCD8Tcell$groupdetail[tumorCD8Tcell$groupdetail%in%c("dMMR-PD1","dMMR-PD1+CTLA4")]<-"dMMR_ICI(+)"
library(Seurat)
library(ggplot2)
library(harmony)
library(SingleCellExperiment)
library(slingshot)
library(scales)
library(viridis)
library(UpSetR)
library(pheatmap)
library(msigdbr)
library(fgsea)
library(knitr)
library(gridExtra)
library(tradeSeq)
library(RColorBrewer)

table(tumorCD8Tcell$Celltype)
sce.HEN.CD4_all <- as.SingleCellExperiment(tumorCD8Tcell, assay = "RNA")

sce.HEN.CD4_all_slinged <- slingshot(sce.HEN.CD4_all, 
                                     clusterLabels = sce.HEN.CD4_all$Celltype, 
                                     reducedDim = reducedDims(sce.HEN.CD4_all)$UMAP, 
                                     start.clus = "CD8+ Tn",
                                     approx_points = 0)
sds <- as.SlingshotDataSet(sce.HEN.CD4_all_slinged)





plot(reducedDims(sce.HEN.CD4_all)$UMAP, col = sce.HEN.CD4_all$color, asp = 1, pch = 16)
lines(as.SlingshotDataSet(sce.HEN.CD4_all_slinged), lwd = 8, col = 'black')
title("CD8+ T cell Lineages in Tumor")

tumorCD8Tcell_CellType_Trajectory_overlay <- DimPlot(tumorCD8Tcell, label = TRUE, label.size = 6, pt.size = 0.75, reduction = "umap", group.by = "Celltype", repel = TRUE, cols = c("#91ccc0", "#7Fabd1", "#F7AC53", "#EC6E66", 
                     "#B5CE4E", "#BD7795")) +
  geom_path(data = as.data.frame(slingCurves(sds)[[1]]$s[slingCurves(sds)[[1]]$ord, ]), aes(umap_1, umap_2), size = 0.5, col = "black") +
  geom_path(data = as.data.frame(slingCurves(sds)[[2]]$s[slingCurves(sds)[[2]]$ord, ]), aes(umap_1, umap_2), size = 0.5, col = "black") +
  ggtitle("CD8+ T cells (subtypes)") + 
  theme(title = element_text(size =15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15))



data <- as.data.frame(slingCurves(sds)[[1]]$s[slingCurves(sds)[[1]]$ord, ])



tumorCD8Tcell_CellType_Trajectory_overlay_NoLabel <- DimPlot(tumorCD8Tcell, label = FALSE, label.size = 6, pt.size = 0.75, reduction = "umap", group.by = "Celltype", repel = TRUE, cols = c("#91ccc0", "#7Fabd1", "#F7AC53", "#EC6E66", 
                                                                                                                                                                                             "#B5CE4E", "#BD7795")) +
  geom_path(data = as.data.frame(slingCurves(sds)[[1]]$s[slingCurves(sds)[[1]]$ord, ]), aes(umap_1, umap_2), size = 1, col = "black") +
  geom_path(data = as.data.frame(slingCurves(sds)[[2]]$s[slingCurves(sds)[[2]]$ord, ]), aes(umap_1, umap_2), size = 1, col = "black")  +
  ggtitle("") + 
  theme(
    text = element_text(family = "Arial"),
    title = element_text(size =15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15))

tumorCD8Tcell_CellType_Trajectory_overlay_NoLabel


ggsave("~/R/ICBstudy/result/Figure/Figure 2/2OandP/CD8_Tcell in tumor.tiff",
       tumorCD8Tcell_CellType_Trajectory_overlay_NoLabel,
       width = 20,       
       height = 15,
       units = "cm",
       dpi = 1800,
       compression = "lzw",
       bg = "white"
)







dMMR <- tumorCD8Tcell@meta.data$groupdetail
table(dMMR)
dMMR[dMMR!= "pMMR-naive"] <- TRUE
table(dMMR)
dMMR[ dMMR == "pMMR-naive"] <- FALSE
table(dMMR)
dMMR <- dMMR == TRUE



sce.HEN.CD4_all_slinged.OnTreatment <- sce.HEN.CD4_all_slinged[, which(dMMR) ]



as.SlingshotDataSet(sce.HEN.CD4_all_slinged.OnTreatment)
as.SlingshotDataSet(sce.HEN.CD4_all_slinged)


####KS-test for the 3 lineages


###Lineage1: CD8+ Tn  CD8+ Tem-CRTAM  CD8+ Tem-GIMAP7  CD8+ Tex  

ks.test(slingPseudotime(sce.HEN.CD4_all_slinged.OnTreatment)[colData(sce.HEN.CD4_all_slinged.OnTreatment)$groupdetail == "dMMR_ICI(-)", 1],
        slingPseudotime(sce.HEN.CD4_all_slinged.OnTreatment)[colData(sce.HEN.CD4_all_slinged.OnTreatment)$groupdetail == "dMMR_ICI(+)", 1], exact = TRUE)


table(colData(sce.HEN.CD4_all_slinged.OnTreatment)$groupdetail)




PT <- slingPseudotime(sce.HEN.CD4_all_slinged.OnTreatment)
#PT[is.na(PT)] <- 0
head(PT)



####Collect the barcodes related to the conditions of interest: Mono-therapy Vs. Combination-therapy
dMMRnaive_Cells<-rownames(tumorCD8Tcell@meta.data[ c(tumorCD8Tcell$groupdetail=="dMMR_ICI(-)"), ])
dMMRPD1_Cells <- rownames(tumorCD8Tcell@meta.data[ c(tumorCD8Tcell$groupdetail=="dMMR_ICI(+)"), ])




length(dMMRnaive_Cells)
length(dMMRPD1_Cells)








###Filter and split by lineage (remove NA values)
PT_1 <- PT[!(is.na(PT[, 1])), 1]





Density_colors <- c("#F7AC53","#91ccc0","#7Fabd1")

####Lineage1: CD8+ Tn  CD8+ Tem-CRTAM  CD8+ Tem-GIMAP7  CD8+ Tex 

ds <- list(dMMRnaive = density(PT_1[names(PT_1) %in% dMMRnaive_Cells]),
           dMMRPD1 = density(PT_1[names(PT_1) %in% dMMRPD1_Cells]))


# Scale the pseudotime from 0 to 100
ds$dMMRnaive$x <- rescale(ds$dMMRnaive$x, to = c(0, 100))
ds$dMMRPD1$x <- rescale(ds$dMMRPD1$x, to = c(0, 100))

xlim <- range(c(ds$dMMRnaive$x, ds$dMMRPD1$x))
ylim <- range(c(ds$dMMRnaive$y, ds$dMMRPD1$y))
###naive vs PD1
plot(xlim, ylim, col = "white", xlab = NULL, ylab = NULL, cex.lab=1, cex.axis=1,yaxt = "n")
axis(2, at = seq(0, 0.12, by = 0.03))
polygon(c(min(ds$dMMRnaive$x), ds$dMMRnaive$x, max(ds$dMMRnaive$x)), c(0, ds$dMMRnaive$y, 0),
  col = alpha(Density_colors[1], alpha = .5))


polygon(c(min(ds$dMMRPD1$x), ds$dMMRnaive$x, max(ds$dMMRnaive$x)), c(0, ds$dMMRPD1$y, 0), # 2 midlle ones wrong?
        col = alpha(Density_colors[2], alpha = .5))





legend("topleft",
       legend = c("dMMR_ICI(-)", "dMMR_ICI(+)"),
       fill = Density_colors[1:2],
       bty = "n", cex = 1, x.intersp = 1, y.intersp = 0.3,
       inset = c(0, -0.14))  








####celltype

Density_colors<-c("#EC6E66","#B5CE4E","#BD7795",  "#52AADC")
table(tumorCD8Tcell$Celltype)
CRTAM<-rownames(tumorCD8Tcell@meta.data[ c(tumorCD8Tcell$Celltype=="CD8+ Tem-CRTAM"), ])
GIMAP7<- rownames(tumorCD8Tcell@meta.data[ c(tumorCD8Tcell$Celltype=="CD8+ Tem-GIMAP7"), ])
Tex<- rownames(tumorCD8Tcell@meta.data[ c(tumorCD8Tcell$Celltype=="CD8+ Tex"), ])
Tn<- rownames(tumorCD8Tcell@meta.data[ c(tumorCD8Tcell$Celltype=="CD8+ Tn"), ])
ds <- list(CRTAM = density(PT_1[names(PT_1) %in% CRTAM]),
           GIMAP7 = density(PT_1[names(PT_1) %in% GIMAP7]),
           Tex = density(PT_1[names(PT_1) %in% Tex]),
           Tn = density(PT_1[names(PT_1) %in% Tn])
)

ds$CRTAM$x <- rescale(ds$CRTAM$x, to = c(0, 100))
ds$GIMAP7$x <- rescale(ds$GIMAP7$x, to = c(0, 100))
ds$Tex$x <- rescale(ds$Tex$x, to = c(0, 100))
ds$Tn$x <- rescale(ds$Tn$x, to = c(0, 100))

xlim <- range(c(ds$CRTAM$x, ds$GIMAP7$x,ds$Tex$x,ds$Tn$x))
ylim <- range(c(ds$CRTAM$y, ds$GIMAP7$y,ds$Tex$y,ds$Tn$y))


plot(xlim, ylim, col = "white", xlab = NULL, ylab = NULL, cex.lab=1, cex.axis=1,yaxt = "n")
axis(2, at = seq(0, 0.5, by = 0.1))



polygon(c(min(ds$Tn$x), ds$Tn$x, max(ds$Tn$x)), c(0, ds$Tn$y, 0),
        col = alpha(Density_colors[1], alpha = .5))


polygon(c(min(ds$CRTAM$x), ds$Tn$x, max(ds$Tn$x)), c(0, ds$CRTAM$y, 0), # 2 midlle ones wrong?
        col = alpha(Density_colors[2], alpha = .5))

polygon(c(min(ds$GIMAP7$x), ds$Tn$x, max(ds$Tn$x)), c(0, ds$GIMAP7$y, 0), # 2 midlle ones wrong?
        col = alpha(Density_colors[3], alpha = .5))
polygon(c(min(ds$Tex$x), ds$Tn$x, max(ds$Tn$x)), c(0, ds$Tex$y, 0), # 2 midlle ones wrong?
        col = alpha(Density_colors[4], alpha = .5))

table(tumorCD8Tcell$Celltype)


legend("topleft",
       legend = c("CD8+ Tn", "CD8+ Tem-CRTAM", "CD8+ Tem-GIMAP7", "CD8+ Tex"),
       fill = Density_colors[1:4],
       bty = "n", cex = 1, x.intersp = 1, y.intersp = 0.3,
       inset = c(0, -0.2))  

