library(Seurat)
library(ggplot2)
library(ggpubr)
tumorCD8Tcell <- readRDS("~/R/ICBstudy/result/data/Tcell/tumorCD8Tcell.rds")
data<-tumorCD8Tcell
data<-data[,data$Celltype!="CD8+ Tn"]
gene1=read.table("progenitor exhausted.txt",header=F)
genes1=as.list(gene1)
data=AddModuleScore(data,features=genes1,name='Progenitor_exhausted')
data$'Progenitor exhausted'<-data$Progenitor_exhausted1





library(Seurat)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)


plot_data <- data@meta.data[, c("Celltype", "Progenitor_exhausted1")]
celltypes <- c("CD8+ Tem-CRTAM", "CD8+ Tem-GIMAP7", "CD8+ Temra", "CD8+ Tex")
plot_data <- plot_data %>% 
  filter(Celltype %in% celltypes) %>%
  mutate(Celltype = factor(Celltype, levels = celltypes))


stat.test <- plot_data %>%
  wilcox_test(Progenitor_exhausted1 ~ Celltype, 
              ref.group = "CD8+ Tem-CRTAM") %>%
  add_xy_position(x = "Celltype", step.increase = 0.20)





stat.test$p.signif <- symnum(stat.test$p, corr = FALSE,
                             cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                             symbols = c("***", "**", "*", "ns"))


p <- ggplot(plot_data, aes(x = Celltype, y = Progenitor_exhausted1)) +
  geom_violin(aes(fill = Celltype), scale = "width", width = 0.8, 
              trim = TRUE, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  

  stat_pvalue_manual(
    stat.test, 
    label = "p.signif", 
    tip.length = 0.01,
    vjust = 0.1,
    size = 9,
    bracket.size = 1.2
  ) +
  
  scale_fill_manual(values = c("#7Fabd1","#F7AC53","#EC6E66","#B5CE4E")) +
  labs(x = NULL, y = "Progenitor Exhausted Score") +
  theme_classic(base_size = 18) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 30),
    axis.text.y = element_text(color = "black", size = 30),
    axis.title.y = element_text(
      size = 30,
      color = "#2c3e50",
      margin = margin(r = 20)
    ),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  ) +
  coord_cartesian(ylim = c(min(plot_data$Progenitor_exhausted1), 
                           max(plot_data$Progenitor_exhausted1) * 1.4))


print(p)
ggsave("Progenitor Exhausted Score.tiff",
       plot = p,
       width = 35,        
       height = 23,
       units = "cm",
       dpi = 1800,
       compression = "lzw",
       bg = "white")
