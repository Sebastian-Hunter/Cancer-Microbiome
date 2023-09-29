library(tidyverse)
library(data.table)
library(ggpubr)
library(vcd)
library("survival")
library("survminer")
setwd("D:/Canada/UBC/UBC Research/Zach_cancer_microbiome")
meta <- read.csv("Metadata-TCGA-All-18116-Samples.csv", header =  TRUE)
reads <- read.csv("Kraken-TCGA-Voom-SNM-Full-Data.csv",header = TRUE,row.names = 1)

breast_cancer <- subset(meta, primary_site == "Breast")
subjects <- breast_cancer$X
interest_bugs <- c("k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Bacillaceae.g__Bacillus", 
                   "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Paenibacillaceae.g__Paenibacillus",
                   "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Staphylococcaceae.g__Staphylococcus")
subject_reads <- reads[subjects,interest_bugs]
table_bugs <- setDT(subject_reads, keep.rownames = TRUE)[]
colnames_bugs <- c("ID", "Bacillus", "Paenibacillus", "Staphylococcus")
colnames(table_bugs) <- colnames_bugs
bugs <- gather(data=table_bugs, genus, value, -ID)
meta2 <- breast_cancer[, c(1,4,6,7,10,17,24,25,26,42)]
meta3 <- meta2[meta2$case_uuid %in% meta2$case_uuid[duplicated(meta2$case_uuid)],]
#filter samples that have normal tissue
meta3 <- meta3%>%group_by(case_uuid)%>%filter(any(sample_type=="Blood Derived Normal" | 
                                                    sample_type=="Solid Tissue Normal"))
#Join variables to simply data
meta3 <- meta3  %>% mutate(stage_label = case_when(
  pathologic_stage_label=="Stage I" | pathologic_stage_label=="Stage IA" | 
    pathologic_stage_label=="Stage IB" ~ "Stage I", pathologic_stage_label=="Stage II" | 
    pathologic_stage_label=="Stage IIA" | pathologic_stage_label=="Stage IIB" ~ "Stage II", 
  pathologic_stage_label=="Stage III" | pathologic_stage_label=="Stage IIIA" | pathologic_stage_label=="Stage IIIB" |
    pathologic_stage_label=="Stage IIIC" ~ "Stage III", pathologic_stage_label=="Stage IV" ~ "Stage IV"),
  t_label = case_when(pathologic_t_label=="T1" | pathologic_t_label=="T1a" | pathologic_t_label=="T1b" |
                        pathologic_t_label=="T1c" ~ "T1", pathologic_t_label=="T2" | pathologic_t_label=="T2a" |
                        pathologic_t_label=="T2b" ~ "T2", pathologic_t_label=="T3" | pathologic_t_label=="T3a" ~ "T3",
                      pathologic_t_label=="T4" | pathologic_t_label=="T4b" | pathologic_t_label=="T4d" ~ "T4"),
  n_label = case_when(pathologic_n_label=="N0" | pathologic_n_label=="N0 (i-)" | pathologic_n_label=="N0 (i+)" | 
                        pathologic_n_label=="N0 (mol+)" ~ "N0", pathologic_n_label=="N1" | pathologic_n_label=="N1a" |
                        pathologic_n_label=="N1b" | pathologic_n_label=="N1c" | pathologic_n_label=="N1mi" ~ "N1",
                      pathologic_n_label=="N2" | pathologic_n_label=="N2a" ~ "N2", pathologic_n_label=="N3" |
                        pathologic_n_label=="N3a" | pathologic_n_label=="N3b" | pathologic_n_label=="N3c" ~ "N3"))  
colnames(meta3)[1] <- "ID"

#merge the whole datset
total <- merge(bugs, meta3, by="ID")

#AGE OF DIAGNOSIS
#prune na values from age_at_diagnosis
total_diag_age <- drop_na(total, age_at_diagnosis)

##Bacillus
ggscatter(subset(total_diag_age, genus =="Bacillus" & sample_type!="Metastatic"), x="value", 
          y="age_at_diagnosis",add="reg.line", conf.int = TRUE, 
          xlab = "Abundance", ylab = "Age at diagnosis", title = "Bacillus") +
  facet_grid(~sample_type)+ stat_cor(
    method = "pearson",
    label.x = 19, label.y = c(80,85))

##Paenibacillus
ggscatter(subset(total_diag_age, genus =="Paenibacillus" & sample_type!="Metastatic"), x="value", 
          y="age_at_diagnosis",add="reg.line", conf.int = TRUE, 
          xlab = "Abundance", ylab = "Age at diagnosis", title = "Paenibacillus") +
  facet_grid(~sample_type)+ stat_cor(
    method = "pearson",
    label.x = 18, label.y = c(80,85,90))

##Staphylococcus
ggscatter(subset(total_diag_age, genus =="Staphylococcus" & sample_type!="Metastatic"), x="value", 
          y="age_at_diagnosis",add="reg.line", conf.int = TRUE, 
          xlab = "Abundance", ylab = "Age at diagnosis", title = "Staphylococcus") +
  facet_grid(~sample_type)+ stat_cor(
    method = "pearson",
    label.x = 17, label.y = c(85,90))

#DAYS TO DEATH
#prune na values from days_to_death
total_days_death <- drop_na(total, days_to_death)

##Bacillus
ggscatter(subset(total_days_death, genus =="Bacillus" & sample_type!="Metastatic"), x="value", 
          y="days_to_death",add="reg.line", conf.int = TRUE, 
          xlab = "Abundance", ylab = "Days to death", title = "Bacillus") +
  facet_grid(~sample_type)+ stat_cor(
    method = "pearson",
    label.x = 16.8, label.y = c(3500,4000))

##Paenibacillus
ggscatter(subset(total_days_death, genus =="Paenibacillus" & sample_type!="Metastatic" & value>14), x="value", 
          y="days_to_death",add="reg.line", conf.int = TRUE, 
          xlab = "Abundance", ylab = "Days to death", title = "Paenibacillus") +
  facet_grid(~sample_type)+ stat_cor(
    method = "pearson",
    label.x = 14, label.y = c(4600,5000))
##Staphylococcus
ggscatter(subset(total_days_death, genus =="Staphylococcus" & sample_type!="Metastatic"), x="value", 
          y="days_to_death",add="reg.line", conf.int = TRUE, 
          xlab = "Abundance", ylab = "Days to death", title = "Staphylococcus") +
  facet_grid(~sample_type)+ stat_cor(
    method = "pearson",
    label.x = 16, label.y = c(4500,5000))

#prune na values from pathologic_t_label
total_t_label <- drop_na(total, t_label)
my_comparisons <- list( c("Blood Derived Normal", "Primary Tumor"), 
                        c("Primary Tumor", "Solid Tissue Normal"), 
                        c("Blood Derived Normal", "Solid Tissue Normal") )
ggplot(subset(total_t_label, genus=="Bacillus" & sample_type!="Metastatic" ), 
       aes(x=sample_type, y=value, color=sample_type))  + geom_boxplot(outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme_bw()+ geom_jitter(position = position_jitterdodge(0.15),alpha=0.4) + facet_grid(.~t_label) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## other comparisons
my_comparisons2 <- list(c("T1", "T2"), c("T2","T3"), c("T3","T4"),c("T1","T3"),c("T2","T4"),c("T1","T4"))
ggplot(subset(total_t_label, genus=="Bacillus" & sample_type!="Metastatic"), 
       aes(x=t_label, y=value, color=sample_type)) + geom_boxplot(outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons2, label="p.signif")+
  theme_bw() + geom_jitter(width = 0.15, alpha=0.4) + facet_wrap(.~sample_type)


#prune na values from pathologic_n_label
total_n_label <- drop_na(total, n_label)
ggplot(subset(total_n_label, genus=="Bacillus" & sample_type!="Metastatic"), 
       aes(x=sample_type, y=value, color=sample_type)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(0.15), alpha=0.4)+
  stat_compare_means(comparisons = my_comparisons, label="p.signif") +
  theme_bw() + facet_grid(.~n_label) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

## other comparisons
my_comparisons3 <- list(c("N0", "N1"), c("N1","N2"), c("N2","N3"),c("N0","N2"),c("N1","N3"),c("N0","N3"))
ggplot(subset(total_n_label, genus=="Bacillus" & sample_type!="Metastatic"), 
       aes(x=n_label, y=value, color=sample_type)) + geom_boxplot(outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons3, label="p.signif")+
  theme_bw()+geom_jitter(position = position_jitterdodge(0.3),alpha=0.4) + facet_wrap(.~sample_type)

#prune na values from pathologic_stage_label
total_stage_label <- drop_na(total, stage_label)
ggplot(subset(total_stage_label, genus=="Bacillus" & sample_type!="Metastatic"), 
       aes(x=sample_type, y=value, color=sample_type)) + geom_boxplot(outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons, label="p.signif")+
  theme_bw()+geom_jitter(position = position_jitterdodge(0.15),alpha=0.4) + facet_grid(.~stage_label)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## other comparisons
my_comparisons4 <- list(c("Stage I", "Stage II"), c("Stage II","Stage III"), c("Stage III","Stage IV"),
                        c("Stage I","Stage III"),c("Stage II","Stage IV"),c("Stage I","Stage IV"))
ggplot(subset(total_stage_label, genus=="Bacillus" & sample_type!="Metastatic" & sample_type!="Blood Derived Normal"), 
       aes(x=stage_label, y=value, color=sample_type)) + geom_boxplot(outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons4, label="p.signif")+
  theme_bw()+geom_jitter(position = position_jitterdodge(0.3),alpha=0.4) + facet_wrap(.~sample_type)

my_comparisons5 <- list(c("Stage I", "Stage II"), c("Stage II","Stage III"),c("Stage I","Stage III"))
ggplot(subset(total_stage_label, genus=="Bacillus" & sample_type=="Blood Derived Normal"), 
       aes(x=stage_label, y=value, color=sample_type)) + geom_boxplot(outlier.shape = NA) + 
  stat_compare_means(comparisons = my_comparisons5,)+
  theme_bw()+geom_jitter(position = position_jitterdodge(0.3),alpha=0.4) 

##Abundance according all cancer markers
total2 <- total %>% group_by(genus, experimental_strategy, n_label, 
                             t_label, stage_label) %>% summarize(value=mean(value))


ggballoonplot(drop_na(subset(total2, experimental_strategy=="RNA-Seq"), 
                      c(n_label,t_label,stage_label)), 
              x="stage_label", 
              y="n_label", fill="value", size="value",
              facet.by = c("genus","t_label"), ggtheme = theme_bw()) + 
  scale_fill_viridis_c(option = "C") + rotate_x_text(angle = 90) 


