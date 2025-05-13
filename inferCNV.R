library(infercnv)
library(Seurat)
colonepi <- readRDS("~/R/ICBstudy/Figure/Figure 1/Infercnv/colonepi.rds")
setwd("~/R/ICBstudy/Figure/Figure 1/Infercnv/result")
outdir_1<-"~/R/ICBstudy/Figure/Figure 1/Infercnv/result"
matrix<-as.matrix(colonepi@assays$RNA@counts)
a<-as.data.frame(colonepi$Type)
a<-as.matrix(a)
colnames(a)<-"cell_type"
gene.order<-read.table(file="gene_position_updated_hg19.txt",header = F,sep = "\t")
rownames(gene.order)<-gene.order$V1
gene.order<-gene.order[,-1]
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix,
                                     annotations_file = a,
                                     gene_order_file = gene.order,
                                     ref_group_names = c("P1-Normal","P2-Normal","P3-Normal",
                                                         "P4-Normal","P5-Normal","P6-Normal",
                                                         "P7-Normal","P8-Normal","P9-Normal"))
infercnvobj <- infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir=outdir_1, 
                              window_length = 101, 
                              max_centered_threshold = 3,
                              plot_steps = T,denoise = TRUE,
                              cluster_by_groups=T,
                              no_prelim_plot = T,
                              HMM=FALSE,
                              sd_amplifier = 1.3,
                              analysis_mode = "samples",
                              num_threads=16,
                              png_res = 2000,
                             write_expr_matrix=T)
save.image(file="infercnv.RData")


infercnv_normal<-read.table("infercnv.references.txt",header = T)
infercnv_malig<-read.table("infercnv.observations.txt",header = T)
a<-as.data.frame(rownames(infercnv_malig),rownames(infercnv_normal))
allinfercnv<- merge(infercnv_malig, infercnv_normal, by = "row.names")
rownames(allinfercnv)<-allinfercnv$Row.names
allinfercnv<-allinfercnv[,-1]
a<-allinfercnv-1
b<-a^2
CNV_score=as.data.frame(colSums(b))
colnames(CNV_score)="CNV_score"
rownames(CNV_score) <- gsub("\\.", "-", rownames(CNV_score))
colonepi$CNVscore <- CNV_score$CNV_score[match(rownames(colonepi@meta.data), rownames(CNV_score))]
data<-colonepi
data$patient<-data$orig.ident
data$patient[data$patient%in%c("C1")]<-"P6"
data$patient[data$patient%in%c("C2")]<-"P7"
data$patient[data$patient%in%c("C3")]<-"P8"
data$patient[data$patient%in%c("C4")]<-"P1"
data$patient[data$patient%in%c("C5")]<-"P3"
data$patient[data$patient%in%c("C6")]<-"P2"
data$patient[data$patient%in%c("C7")]<-"P4"
data$patient[data$patient%in%c("C8")]<-"P5"
data$patient[data$patient%in%c("C9")]<-"P9"
table(data$groupdetail,data$patient)
data$type1<-data$original
data@meta.data<- within(data@meta.data, {
  type1<- sapply(strsplit(original, "_"), `[`, 2)
})

data$sample <- paste(data$patient, data$type1, sep = "_")

table(data$sample)
library(ggplot2)
table(data$patient)

table(data$type1,data$Type)
normaldata<-data[,data$type1=="T"]
table(normaldata$groupdetail)
normaldata$group<-normaldata$groupdetail
normaldata$group[normaldata$groupdetail%in%c("dMMR-PD1", "dMMR-PD1+CTLA4")]<-"dMMR_ICI(+)"
normaldata$group[normaldata$groupdetail%in%c("dMMR-naive")]<-"dMMR_ICI(-)"
normaldata$group[normaldata$groupdetail%in%c("pMMR-naive")]<-"pMMR_ICI(-)"
table(normaldata$group)
normaldata$group<-factor(normaldata$group,levels = c("pMMR_ICI(-)","dMMR_ICI(-)", "dMMR_ICI(+)"))



library(ggplot2)
library(ggpubr)


normaldata$group <- factor(normaldata$group, levels = c("pMMR_ICI(-)","dMMR_ICI(-)", "dMMR_ICI(+)"))
table(normaldata$group)
ggplot(normaldata@meta.data, aes(x = group, y = CNVscore, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") + 
  stat_compare_means(comparisons = list(c("pMMR_ICI(-)","dMMR_ICI(-)"),
                                        c("pMMR_ICI(-)","dMMR_ICI(+)"),
                                        c("dMMR_ICI(-)", "dMMR_ICI(+)")), 
                     method = "wilcox.test", label = "p.signif", size = 9,bracket.size = 1.2,step.increase =0.5 ) +  
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73")) + 
  labs(x = NULL, y = "CNV Score") + 
  theme_classic() +  
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),  
    legend.position = "none",  
    plot.title = element_text(size = 18, hjust = 0.5) 
  )+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20))


ggsave("CNV score.tiff",
       width = 15,        
       height = 18,
       units = "cm",
       dpi = 3600,
       compression = "lzw",
       bg = "white")  


ggplot(normaldata@meta.data, aes(x = group, y = CNVscore, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, color = "black") +  
  stat_compare_means(comparisons = list(c("pMMR_ICI(-)","dMMR_ICI(-)"),
                                        c("pMMR_ICI(-)","dMMR_ICI(+)"),
                                        c("dMMR_ICI(-)", "dMMR_ICI(+)")), 
                     method = "wilcox.test", label = "p.signif", size = 9,bracket.size = 1.2,step.increase =0.5 ) +  
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73")) + 
  labs(x = NULL, y = "CNV Score") +  
  theme_classic() +  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 14, color = "black"),   
    axis.title = element_text(size = 16),    
    legend.position = "none",  
    plot.title = element_text(size = 18, hjust = 0.5) 
  )+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100))


ggsave("CNV score1.tiff",
       width = 15,         
       height = 18,
       units = "cm",
       dpi = 3600,
       compression = "lzw",
       bg = "white")  



