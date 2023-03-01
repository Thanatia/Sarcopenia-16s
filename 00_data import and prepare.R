if(T){library(phyloseq)
  library(readxl)
  library(stringr)
  library(dplyr)
  library(phyloseq)
  library(vegan)
  library(Biostrings) # To read fasta file
  library(reshape2)
  library(ape) # to read tree file
  library(scales)
  library(picante)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
  library(ggbiplot)
  library(pairwiseAdonis)
  library(ggtree)
  library(ape)
  library(biomformat)
  library(reshape2)
  library(RColorBrewer)
  library(lars)
  library(glmnet)
  library(pROC) 
  library(randomForest)
}
Sarc_tree <- read.tree('./phylogeny_tree/phylogeny.tre')
ASV_Data <- read.csv(file = './ASV.csv')
rownames(ASV_Data) <- ASV_Data$OTU.ID
ASV_Data <- ASV_Data[,-1]
PD_Diag <- read.csv(file = './Huaxi_diag.csv')
PD_Mapping <- read.csv(file = './ID_mapping.csv')
PD_Mapping <- PD_Mapping[,1:5]
PD_Blood <- read.csv(file = './Blood_cond.csv')
PD_Meta <- read.csv(file = './Patient_meta.csv')
PD_merge <- inner_join(PD_Mapping,PD_Diag,by = 'Huaxi_ID')
PD_Clin <- inner_join(PD_Meta,PD_merge,by = 'Huaxi_ID')
rownames(PD_Clin) <- PD_Clin$My_ID
data.frame(seq(1:ncol(PD_Clin)),colnames(PD_Clin),t(PD_Clin[1,]))
PD_Clin <- PD_Clin[order(PD_Clin$MY_ID),]
head(PD_Clin) 
PD_Clin_fixed <- PD_Clin[,-c(3,5,50,52,54,55,56,57)]
PD_Clin_fixed$Sarcopenia <- ifelse(PD_Clin_fixed$Diag==0,'Normal','Sarcopenia')
BMI <- as.numeric(PD_Clin_fixed$weight)/((as.numeric(PD_Clin_fixed$height)/100)^2)
hist(BMI)
PD_Clin_fixed$BMI <- BMI
saveRDS(PD_Clin_fixed,file = 'Sarc_pd.rds')
PD_Clin_fixed <- readRDS('Sarc_pd.rds')
saveRDS(PD_Blood,file = 'Sarc_blood.rds')
PD_Blood <- readRDS('Sarc_blood.rds')
# taxonomy====
dat.b<-read_biom("./01.ASV.tax/02.normalize/asv_taxa_table.biom")
mp0 = import_biom(dat.b, parseFunction = parse_taxonomy_greengenes)
tax_table(mp0) <- tax_table(mp0)[, 1:7]
sample_names(mp0) <- colnames(ASV_Data)
ASV_tax <- tax_table(mp0)
ASV_meta <- sample_data(PD_Clin_fixed)
Sarc_ASV = merge_phyloseq(mp0, ASV_meta)
saveRDS(Sarc_ASV,file = 'Sarc_ASV.rds')
Sarc_ASV <- readRDS('Sarc_ASV.rds')
Diag <- PD_Clin_fixed$Sarcopenia
Sarc_ASV <- merge_phyloseq(Sarc_ASV,Sarc_tree)
# ASV filter ====
total = median(sample_sums(Sarc_ASV))
standf = function(x, t=total) round(t * (x / sum(x)))
gps = transform_sample_counts(Sarc_ASV, standf)
sample_sums(gps)
dim(otu_table(gps))
gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)
gpsf_per <- transform_sample_counts(gps, function(x) x / sum(x) )
dim(otu_table(gpsf))
GP = filter_taxa(gpsf, function(x) sum(x > 2) > (0.05*length(x)), TRUE)
dim(otu_table(GP))
GPr  = transform_sample_counts(GP, function(x) x / sum(x) )
GPr 
otu_table(GPr)[1:5][1:5]
GPfr = filter_taxa(gpsf_per, function(x) mean(x) > 1e-5, TRUE)
dim(otu_table(GPfr))
Sarc_ASV_fr <- GPfr
saveRDS(Sarc_ASV_fr,file = 'Sarc_ASV_fr.rds')
SARC_ASV_counts <- otu_table(Sarc_ASV)
Normal_ASV_counts <- SARC_ASV_counts[,Sarc_diag=='Normal']
Sarc_ASV_counts <- SARC_ASV_counts[,Sarc_diag!='Normal']
shuling <- function(x){sum(x>0)}
dim(Sarc_ASV_counts)
Sarc_ASV_id <- rownames(Sarc_ASV_counts)[apply(Sarc_ASV_counts, 1, shuling)>1]
Norm_ASV_id <- rownames(Normal_ASV_counts)[apply(Normal_ASV_counts, 1, shuling)>1]
library(venn)
library(venn)
a1<-list(Sarc_ASV_id, Norm_ASV_id)
names(a1) <- c('Sarcopenia','Normal')
venn(a1,zcolor='style')
pdf(file="venn.pdf",width=8,height=8)