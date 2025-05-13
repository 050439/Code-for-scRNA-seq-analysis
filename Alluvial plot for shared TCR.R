rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(tidyr)
tumorCD8Tcell <- readRDS("~/R/ICBstudy/result/data/Tcell/tumorCD8Tcell.rds")
LNCD8T <- readRDS("~/R/ICBstudy/result/data/Tcell/LNCD8T.rds")

table(LNCD8T$Celltype)
LN<-LNCD8T[,!LNCD8T$Celltype%in%c("MAIT","CD8+ Tn")]
tumor<-tumorCD8Tcell[,!tumorCD8Tcell$Celltype%in%c("MAIT","CD8+ Tn")]
table(LN$TCR)
LN<-LN[,LN$TCR%in%c("A","AAB","AB","B")]
tumor<-tumor[,tumor$TCR%in%c("A","AAB","AB","B")]

LN<-LN[,LN$groupdetail!="pMMR-naive"]
LN$groupdetail[LN$groupdetail=="dMMR-naive"]<-"dMMR_ICI(-)"
LN$groupdetail[LN$groupdetail%in%c("dMMR-PD1","dMMR-PD1+CTLA4")]<-"dMMR_ICI(+)"
tumor<-tumor[,tumor$groupdetail!="pMMR-naive"]
tumor$groupdetail[tumor$groupdetail=="dMMR-naive"]<-"dMMR_ICI(-)"
tumor$groupdetail[tumor$groupdetail%in%c("dMMR-PD1","dMMR-PD1+CTLA4")]<-"dMMR_ICI(+)"




LN_data <- LN@meta.data %>%
  dplyr::select(c(patient,Celltype, CTnt, groupdetail)) %>%
  filter(!is.na(CTnt)) %>%
  mutate(Sample = "LN")
tumor_data <- tumor@meta.data %>%
  dplyr::select(c(patient,Celltype, CTnt, groupdetail)) %>%
  filter(!is.na(CTnt)) %>%
  mutate(Sample = "Tumor")
All<- bind_rows(LN_data, tumor_data)

All<- All %>% rename(Patient = patient, Sampletype = Sample,Treatment=groupdetail)
All<- All %>% rename(cdr3s_nt = CTnt)
All<- All %>% rename(annotation = Celltype)
samples <- All%>% select(c(Patient, Sampletype, Treatment)) %>% distinct()
TCR_table <- All %>% filter(!is.na(cdr3s_nt)) %>% mutate(pt_cdr3nt = paste0(Patient, "_", cdr3s_nt)) %>% tibble::rownames_to_column("seurat_barcode") %>% select(c(seurat_barcode, pt_cdr3nt, Patient, Treatment, Sampletype, annotation)) 
TCR<-TCR_table
TME <- TCR %>% filter(Sampletype == "Tumor")
LN <- TCR %>% filter(Sampletype == "LN")
TME_sh <- TME %>% filter(pt_cdr3nt %in% LN$pt_cdr3nt)
n_distinct(TME_sh$pt_cdr3nt)

LN_sh <- LN %>% filter(pt_cdr3nt %in% TME$pt_cdr3nt)
n_distinct(LN_sh$pt_cdr3nt)


TME_sh_freq <- TME_sh %>%
  group_by(pt_cdr3nt, annotation, Treatment) %>%
  count() %>%
  rename(TME = annotation, Freq_TME = n) %>%
  ungroup()

LN_sh_freq <- LN_sh %>%
  group_by(pt_cdr3nt, annotation,  Treatment) %>%
  count() %>%
  rename(LN = annotation,  Freq_LN = n) %>%
  ungroup()

Tcellsharing <- TME_sh_freq %>% full_join(LN_sh_freq, by = c("pt_cdr3nt", "Treatment"),relationship = "many-to-many")
Tcellsharing3 <- Tcellsharing %>% select(!c(Freq_TME, Freq_LN)) %>% group_by(TME,  LN,  Treatment) %>% summarize(Freq_TME = n_distinct(pt_cdr3nt), Freq_LN = n_distinct(pt_cdr3nt))

TME <- TCR %>% filter(Sampletype == "Tumor")
LN <- TCR %>% filter(Sampletype == "LN")
patients_2samples <- intersect(LN$Treatment, TME$Treatment)
patients_2samples
ann <- LN %>% group_by(Treatment, annotation) %>% 
  summarize(nr_LN = n_distinct(pt_cdr3nt)) 
tot_ln <- LN %>% group_by(Treatment) %>% summarize(tot_LN = n_distinct(pt_cdr3nt)) 
nr_TCR_LN <- ann %>% left_join(tot_ln) %>% dplyr::rename(LN = annotation)

#TME
ann <- TME  %>% group_by(Treatment, annotation) %>% 
  summarize(nr_TME = n_distinct(pt_cdr3nt)) 
tot_TME <- TME  %>% group_by(Treatment) %>% summarize(tot_TME = n_distinct(pt_cdr3nt)) 
nr_TCR_TME <- ann %>% left_join(tot_TME) %>% dplyr::rename(TME = annotation)

Tcellsharing3$combination <- paste(Tcellsharing3$LN, Tcellsharing3$TME, sep = "_")

tmp <- Tcellsharing3 %>% full_join(tot_ln) %>% full_join(tot_TME) %>% ungroup


patient_combination <- tmp %>% tidyr::expand(combination, Treatment) %>% filter(!is.na(combination))
combination_info <- Tcellsharing3 %>% select(LN, TME,  combination) %>% distinct()     
patient_combination <- patient_combination %>% left_join(combination_info)
prop <- Tcellsharing3 %>% full_join(patient_combination) %>% 
  tidyr::replace_na(list(Freq_LN=0, Freq_TME=0)) %>% 
  left_join(tot_ln) %>% left_join(tot_TME) %>% 
  mutate(Prop_LN = Freq_LN / tot_LN, Prop_TME = Freq_TME / tot_TME)      

prop_group <- prop %>% dplyr::group_by(combination, Treatment) %>% 
  mutate(mean_prop_TME = 100*sum(Prop_TME, na.rm = T), mean_prop_LN = 100*sum(Prop_LN, na.rm = T)) %>% ungroup() %>% 
  select(!c(  Freq_TME, Freq_LN, tot_LN, tot_TME, Prop_LN, Prop_TME)) %>% distinct() %>% arrange(combination)

long <- to_lodes_form(prop_group, key = "loc", axes = c("LN", "TME"))


sum(long$mean_prop_LN)
df <- long %>% rename(TME = mean_prop_TME, LN = mean_prop_LN) %>% tidyr::pivot_longer(cols = c(TME, LN)) %>% filter(loc == name) %>% rename(mean_prop = value) %>% select(!name)
table(df$loc)
df$loc<-as.character(df$loc)
df$loc[df$loc=="LN"]<-"TDLN"
df$loc[df$loc=="TME"]<-"Tumor"
df$loc<-factor(df$loc,levels = c("TDLN","Tumor"))

df <- df  %>% mutate(mean_prop_0 = tidyr::replace_na(mean_prop, 0))
sum(df$mean_prop)

colors<-c("#91ccc0","#7Fabd1","#F7AC53","#EC6E66","#B5CE4E","#BD7795","#7C7979","#963B79","#97D0C5","#F39865",
          "#52AADC","#C7C1DE","#EEB6D4","#C89736","#2D8875","#D75B4E")

plot <- ggplot(data = df, aes(x = loc, y = mean_prop_0, stratum = stratum, 
                              alluvium = alluvium, fill = stratum, label = stratum)) +
  geom_flow(stat = "flow", knot.pos = 1/4, aes.flow = "forward") + 
  geom_stratum() +
  scale_fill_manual(values = colors) + 
  xlab("") +
  ylab("Proportion of shared TCR clonotypes (%)") +
  ggtitle("") +
  theme(
    text = element_text(size = 18),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 16),
    legend.title = element_blank(),
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.text = element_text(size = 26, colour = "black"), 
    axis.ticks.x = element_blank(), 
    strip.background = element_blank(), 
    panel.border = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.spacing.y = unit(0.5, "cm")
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  coord_cartesian(ylim = c(0, 70))

print(plot + facet_grid(~Treatment) + theme(legend.position = "right"))

ggsave("TCRsharing between tissue.tiff",
       width = 18,         
       height = 15,
       units = "cm",
       dpi = 3600,
       compression = "lzw",
       bg = "white") 






