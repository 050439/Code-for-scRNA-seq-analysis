library(Seurat)

library(dplyr)

tumorCD8Tcell <- readRDS("~/R/ICBstudy/result/data/Tcell/tumorCD8Tcell.rds")
tumorCD8Tcell<-tumorCD8Tcell[,tumorCD8Tcell$TCR%in%c("A","B","AB","AAB")]
tumorCD8Tcell$TCRspecific<-paste(tumorCD8Tcell$patient,tumorCD8Tcell$CTnt,sep = "_")
LNCD8T <- readRDS("~/R/ICBstudy/result/data/Tcell/LNCD8T.rds")
LNCD8T<-LNCD8T[,LNCD8T$TCR%in%c("A","B","AB","AAB")]
LNCD8T$TCRspecific<-paste(LNCD8T$patient,LNCD8T$CTnt,sep = "_")
share<-intersect(tumorCD8Tcell$TCRspecific,LNCD8T$TCRspecific)

LNCD8T<-LNCD8T[,LNCD8T$Celltype=="CD8+ Tem"]
LNCD8T$share<-LNCD8T$TCRspecific
LNCD8T$share[LNCD8T$TCRspecific%in%share]<-"Share"
LNCD8T$share[!LNCD8T$TCRspecific%in%share]<-"Nonshare"
table(LNCD8T$share)
myeloid <- readRDS("~/R/ICBstudy/result/data/myeloid/myeloid.rds")
myeloidLN<-myeloid[,myeloid$Type=="LN"]
table(Idents(myeloidLN),myeloidLN$groupdetail)
myeloidLN$Celltype<-Idents(myeloidLN)

myeloidLN<-myeloidLN[,myeloidLN$Celltype%in%c("pDC","cDC1","cDC2","DC-LAMP3")]

LNCD8T$Celltype<-as.character(LNCD8T$Celltype)
LNCD8T$Celltype[LNCD8T$share=="Share"]<-"Shared CD8+ Tem"
LNCD8T$Celltype[LNCD8T$share!="Share"]<-"Nonshared CD8+ Tem"
LN<- merge(myeloidLN, y = LNCD8T)  
Idents(LN)<-LN$Celltype


library(stringr)
library(Seurat)
library(CellChat)
library(dplyr)
library(stringr)
control<-LN


control.input <- GetAssayData(control, assay = "RNA", slot = "data") # controlized data matrix
labels <- factor(control$Celltype,levels=levels(Idents(control)))
meta <- data.frame(group = labels, row.names = rownames(control@meta.data)) # create a dataframe of the cell labels
cellchat_control <- createCellChat(object = control.input, meta = meta, group.by = "group")



cellchat_control@DB <- CellChatDB.human
cellchat_control <- subsetData(cellchat_control) # subset the expression data of signaling genes for saving computation cost
cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
cellchat_control <- projectData(cellchat_control, PPI.human)
cellchat_control <- computeCommunProb(cellchat_control, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_control <- filterCommunication(cellchat_control, min.cells = 1)
cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_control <- netAnalysis_computeCentrality(cellchat_control, slot.name = "netP")
cellchat_control <- aggregateNet(cellchat_control)
groupSize <- as.numeric(table(cellchat_control@idents))
groupSize


?netVisual_circle

netVisual_circle(cellchat_control@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "",
                 sources.use=c( "cDC1", "cDC2","DC-LAMP3","pDC"),
                 targets.use =c("Shared CD8+ Tem", "Nonshared CD8+ Tem"))

netVisual_circle(cellchat_control@net$count, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "",
                 sources.use=c( "cDC1", "cDC2","DC-LAMP3","pDC"),
                 targets.use =c("Shared CD8+ Tem", "Nonshared CD8+ Tem"))



netVisual_bubble(cellchat_control, sources.use=c( "cDC1", "cDC2","DC-LAMP3","pDC"),
                 targets.use =c("Shared CD8+ Tem", "Nonshared CD8+ Tem"), remove.isolate = FALSE)









