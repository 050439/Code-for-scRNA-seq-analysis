rm(list = ls())
gc()
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
table(myeloid$Celltype)
myeloidLN<-myeloidLN[,myeloidLN$Celltype=="cDC1"]
LNCD8T$Celltype<-as.character(LNCD8T$Celltype)
LNCD8T$Celltype[LNCD8T$share=="Share"]<-"Shared-CD8+ Tem"
LNCD8T$Celltype[LNCD8T$share!="Share"]<-"Nonshared-CD8+ Tem"
LNCD8T<-LNCD8T[,LNCD8T$Celltype=="Shared-CD8+ Tem"]
LN<- merge(myeloidLN, y = LNCD8T)  
LN<-LN[,LN$groupdetail!="pMMR-naive"]
LN$groupdetail[LN$groupdetail=="dMMR-naive"]<-"-ICI"
LN$groupdetail[LN$groupdetail%in%c("dMMR-PD1","dMMR-PD1+CTLA4")]<-"+ICI"


seuratObj<-LN


library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
seuratObj <- UpdateSeuratObject(seuratObj)
seuratObj <- alias_to_symbol_seurat(seuratObj, "human")
class(seuratObj)


ligand_target_matrix <- readRDS("~/R/ICBstudy/analysedata/niche/human/ligand_target_matrix_nsga2r_final.rds")
lr_network<- readRDS("~/R/ICBstudy/analysedata/niche/human/lr_network_human_21122021.rds")
weighted_networks<- readRDS("~/R/ICBstudy/analysedata/niche/human/weighted_networks_nsga2r_final.rds")
lr_network <- lr_network %>% distinct(from, to)

seuratObj$Celltype<-Idents(seuratObj)
table(seuratObj$Celltype)
seuratObj$Celltype<-as.character(seuratObj$Celltype)
Idents(seuratObj)<-seuratObj$Celltype
table(seuratObj$Celltype)
receiver = "CD8+ Tem"
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

table(seuratObj$Celltype)
sender_celltypes <- c("cDC1")

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

seuratObj$Group<-seuratObj$groupdetail
table(seuratObj$Group)


condition_oi <-  "+ICI"
condition_reference <- "-ICI"
seurat_obj_receiver <- subset(seuratObj, idents = receiver)
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "Group",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]


background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.95) 
order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")




ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))


p_dotplot <- DotPlot(subset(seuratObj, Celltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot



ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()

ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "#EC6E66", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "#EC6E66")

p_ligand_target

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

