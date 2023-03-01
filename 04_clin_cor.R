library(CBCgrps)
library(caret)
p2star <- function(p){
  t1 <- symnum(p,cutpoints = c(0,0.0001,0.001,0.01,0.055,1),symbols = c('****','***','**','*',''),na = NA)
  t1 <- as.character(t1)
  return(t1)
}
# clinical_data ====
clin_comp <- data.frame(clin_data)
clin_comp_sel <- as.data.frame(clin_comp[,-c(1,3,5,9,10,12,13,44,46,48,49,50)])#数值型
tab1 <- twogrps(clin_comp_sel, gvar = "Diag")
clin_3line <- tab1$Table
colnames(clin_3line) <- clin_3line[1,]
clin_3line$p[clin_3line$p == "< 0.001"] <- 0.001
clin_3line$p <- as.numeric(clin_3line$p)
clin_sig <- str_split_fixed(na.omit(clin_3line[clin_3line$p<0.05,])[,1],',',2)[,1]
clin_sig_data <- clin_comp %>% select(clin_sig)
clin_comp_sel <- clin_comp[,-c(1,4)]
tab1 <- twogrps(clin_comp_sel, gvar = "Diag",workspace=2e8)
write.csv(tab1$Table,'clin_sanxian.csv') #临床数据
# blood_data ====
blood_cond <- read.csv('./Blood_cond.csv')
blood_sel <- blood_cond[-112,-c(1,2,4)]
blood_sel$BAS_per <- blood_sel$BAS_per+runif(nrow(blood_sel))/100
blood_sel$BAS_abs <- blood_sel$BAS_abs+runif(nrow(blood_sel))/100
blood_sel <- blood_cond[-112,-c(1,2,4)]
tab2 <- twogrps(blood_sel, gvar = "Diag",workspace = 2e8)
blood_3line <- tab2$Table
colnames(clin_3line) <- clin_3line[1,]
blood_3line$p <- as.numeric(blood_3line$p)
library(stringr)
clin_sig <- str_split_fixed(na.omit(clin_3line[clin_3line$p<0.05,])[,1],',',2)[,1]
blood_sig <- str_split_fixed(na.omit(blood_3line[blood_3line$p<0.05,])[,1],',',2)[,1]
blood_sig_dat <- blood_cond[,blood_sig]
na_num <- which(is.na(blood_sig_dat$VD3))
library(tidyr)
library(tidyverse)
clin_sig_data <- clin_comp %>% select(clin_sig)
blood_sig_data <- blood_sel %>% select(blood_sig)
blood_sig_data_nona <- na.omit(blood_sig_data)
dim(blood_sig_data_nona) 
rownames(blood_sig_data) <- blood_cond$MY_ID[-112]
blood_sig_data_nona <- na.omit(blood_sig_data)
asv_blood_pare <- asv_sel[rownames(blood_sig_data_nona),]
# correlation test====
sel_ASV_box2 <- LEFSE_sig
sel_ASV_box_tocor <- as.data.frame(t(sel_ASV_box2))
clin_sig_data_all_sel <- clin_sig_data_sel[rownames(sel_ASV_box_tocor),]
clin_sig_data_sarc_only <- clin_sig_data_sel[Sarc_diag=='Sarcopenia',]# patients
clin_sig_data_norm_only <- clin_sig_data_sel[Sarc_diag!='Sarcopenia',]# NC
library(corrplot)
color_1<-colorRampPalette(c("#063062","white","red"))#色板
?corrplot
corrplot(cor(clin_sig_data_sarc_only),col = color_1(255), addCoef.col = 'darkblue',tl.col='black',title = 'Sarcopenia',mar = c(2,2,2,2),tl.cex = 1.5,number.cex = 1.5)
corrplot(cor(clin_sig_data_norm_only),col = color_1(255), addCoef.col = 'darkblue',tl.col='black',title = 'Normal',mar = c(2,2,2,2),tl.cex = 1.5,number.cex = 1.5)
asv_sel_sarc_only <- sel_ASV_box_tocor[rownames(clin_sig_data_sarc_only),]
dim(asv_sel_sarc_only)
asv_sel_norm_only <- sel_ASV_box_tocor[rownames(clin_sig_data_norm_only),]
dim(asv_sel_norm_only)
library(psych)
clin_asv_cor <- corr.test(asv_sel_norm_only,clin_sig_data_norm_only,method = 'spearman',adjust = 'none')
clin_asv_cor <- corr.test(asv_sel_sarc_only,clin_sig_data_sarc_only,method = 'spearman',adjust = 'none')
clin_asv_cor_all <- corr.test(sel_ASV_box_tocor,clin_sig_data_all_sel,method = 'spearman',adjust = 'none')
clin_asv_cor_r <- clin_asv_cor$r
clin_asv_cor_p <- clin_asv_cor$p
all_r <- clin_asv_cor_r
all_p <- clin_asv_cor_p
View(all_p)
all_star <- apply(all_p,2,p2star)
rownames(all_star) <- rownames(all_r)
sarc_r2 <- clin_asv_cor_r
sarc_p2 <- clin_asv_cor_p
sarc_psatr2 <- apply(sarc_p2,2,p2star)
norm_r2 <- clin_asv_cor_r
norm_p2 <- clin_asv_cor_p
norm_pstar2 <- apply(norm_p2,2,p2star)
all_r_sel <- all_r[rowMin(all_p)<0.055,]
all_star_sel <- all_star[rownames(all_r_sel),]
pheatmap(norm_r2,display_numbers = norm_pstar2,border_color = NA,cellwidth = 10,cellheight = 10,
         fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = F,cluster_cols = F,main = 'normal')
pheatmap(sarc_r2,display_numbers = sarc_psatr2,border_color = NA,cellwidth = 10,cellheight = 10,
         fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = F,cluster_cols = F,main = 'sarcopenia')
pheatmap(all_r,display_numbers = all_star,border_color = NA,cellwidth = 10,cellheight = 10,
         fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = F,cluster_cols = F,main = 'All_patients')
ph1 <- pheatmap(all_r_sel,display_numbers = all_star_sel,border_color = NA,cellwidth = 25,cellheight = 15, annotation_row = annot_row,annotation_colors = ant_color,
                fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = T,cluster_cols = T,main = 'All_patients')
sel_rown <- str_split_fixed(rownames(all_r_sel)[ph1$tree_row$order],'_',2)[,1]
ann_colors <- list(Phylum = c(Control="#00a0e9",Cachexia='#ff5252'))
mycol_13 <- c(mycol_12,'grey50') 
unique(annot_row$Order)
o1 <- str_c(unique(annot_row$Order),mycol_13,sep = '"="')
unique(annot_row$Phylum)
ant_color = list(
  Order = c(
    "Lachnospirales"="#FAFD7CFF" ,   "Clostridia UCG-014"="#82491EFF","Oscillospirales"="#24325FFF",  
    "Clostridiales"="#B7E4F9FF",     "Christensenellales"="#FB6467FF","Erysipelotrichales"="#526E2DFF",     
    "Lactobacillales"="#E762D7FF",   "Acidaminococcales"="#E89242FF", "Peptostreptococcales-Tissierellales" = "#FAE48BFF",
    "Bifidobacteriales"="#A6EEE6FF", "Coriobacteriales"="#917C5DFF",  "Burkholderiales"="#69C8ECFF",  
    "Bacteroidales"="grey50"  ),
  Phylum = c(
    "Firmicutes" = '#3F51B5', 
    "Actinobacteriota"= '#C62828',
    "Proteobacteria"=  '#009688',
    "Bacteroidota"= '#F44336'   
    
  ) )
all_r_sep<- cbind(norm_r2,sarc_r2)
all_p_sep <- cbind(norm_p2,sarc_p2)
all_r_sep_sel <- all_r_sep[str_split_fixed(rownames(all_r_sel),'_',2)[,1],]
all_p_sep_sel <- all_p_sep[str_split_fixed(rownames(all_r_sel),'_',2)[,1],]
fix(all_p_sep_sel)
table(sarc_p2<0.05)
table(norm_p2<0.05)
all_r_sep_sel[all_p_sep_sel>0.05] <- NA
Nvalue <- apply(all_r_sep_sel, 1, function(x)
{sum(is.na(x))})
all_r_sep_sel <- all_r_sep_sel[!Nvalue==14,]
r_di_plot <- as.data.frame(t(all_r_sep_sel))
r_di_plot$group <- c(rep('Normal',7),rep('Sarcopnenia',7))
r_di_plot$Para <- rownames(r_di_plot)
r_di_plot$Para <- str_split_fixed(r_di_plot$Para,'\\.',2)[,1]
r_ar_di_plotll_plot <- melt(r_di_plot, id = c('group','Para'))
colnames(r_ar_di_plotll_plot)[3] <- 'ASV'
head(r_ar_di_plotll_plot)
r_ar_di_plotll_plot$value <- r_ar_di_plotll_plot$value*1.5
asv_value <- data.frame(
  ASV = rownames(sel_ASV_box2),
  Counts = rowMeans(sel_ASV_box2))
library(dplyr)
r_ar_di_plotll_plot<- inner_join(r_ar_di_plotll_plot,asv_value)
r_ar_di_plotll_plot$Counts <- ifelse(is.na(r_ar_di_plotll_plot$value), NA, r_ar_di_plotll_plot$Counts)
r_ar_di_plotll_plot$Counts <- log2(r_ar_di_plotll_plot$Counts)
unique(r_ar_di_plotll_plot$Para)
asv_le <- intersect(sel_rown,unique(r_ar_di_plotll_plot$ASV))
r_ar_di_plotll_plot$Para <- factor(r_ar_di_plotll_plot$Para,levels = c("t121" ,  "t122" , "lgrip",  "rgrip" , "grip" , "weight", "BMI"))
r_ar_di_plotll_plot$ASV <- factor(r_ar_di_plotll_plot$ASV, levels = rev(asv_le))
ASV_ful <- rownames(all_r_sel)[str_split_fixed(rownames(all_r_sel),'_',2)[,1]%in%asv_le]
ASV_ful <- ASV_ful[match(asv_le,str_split_fixed(ASV_ful,'_',2)[,1])]
ASV_index <- data.frame(
  ASV_tax = ASV_ful,
  ASV = asv_le
)
r_ar_di_plotll_plot2 <- inner_join(r_ar_di_plotll_plot,ASV_index)
r_ar_di_plotll_plot2$ASV_tax <- factor(r_ar_di_plotll_plot2$ASV_tax,levels = rev(ASV_index$ASV_tax))
r_ar_di_plotll_plot2 <- r_ar_di_plotll_plot2[!str_detect(r_ar_di_plotll_plot2$ASV_tax,'uncultured'),]
r_ar_di_plotll_plot2 <- r_ar_di_plotll_plot2[!str_detect(r_ar_di_plotll_plot2$ASV_tax,'norank'),]
r_ar_di_plotll_plot2$value <- r_ar_di_plotll_plot2$value/2
?scale_color_continuous
show_point_shapes()
p1 <- ggplot(r_ar_di_plotll_plot2,aes(x = Para, y = ASV_tax, shape = group, size = Counts))+geom_point(aes(color = value))+theme_bw()+scale_shape_manual(values = c(15,16))
facet(p1 , facet.by = "group",scales="free_y",ncol = 2)+scale_color_gradient2(low = "blue", midpoint = 0, mid = "white", high = "red")+
  # stat_compare_means(aes(group = group),label = "p.format",
  #                    vjust = 1, method = 'wilcox.test', size = 4)+
  # scale_color_brewer(palette = "Set1")+
  theme_classic()
Sarc_asv_rho <- Sarc_ASV_tax[rownames(all_r_sel),f]
Sarc_asv_rho_hasname <- str_c(rownames(Sarc_asv_rho),Sarc_asv_rho$Genus,sep = '_')
Sarc_family_rho <- Sarc_asv_rho$Family
Sarc_order_rho <- Sarc_asv_rho$Order
identical(rownames(all_r_sel), rownames(Sarc_asv_rho))
rownames(all_r_sel) <- Sarc_asv_rho_hasname
annot_row <- data.frame(
  row.names = rownames(all_r_sel),
  Order = Sarc_order_rho,
  Phylum = Sarc_asv_rho$Phylum
)
clin_sig_data_sarc_only
dim(clin_sig_data_sel)
corr.test(clin_sig_data_sel)
clin_sig_data_sel
# asv_sel <- data.frame(otu_table(Sarc_ASV))[,-112]
# asv_sel <- as.data.frame(t(asv_sel))
dim(sel_ASV_box2)
rownames(sel_ASV_box2)
asv_sel <- as.data.frame(t(sel_ASV_box2[,-122]))
dim(asv_sel)
# asv_sel_sarc_only <- asv_sel[rownames(clin_sig_data_sarc_only),]
# asv_sel_norm_only <- asv_sel[rownames(clin_sig_data_norm_only),]
# asv_sel_sarc_only_tocor <- asv_sel_sarc_only[,colMaxs(as.matrix(asv_sel_sarc_only))>100]
# asv_sel_sarc_only_tocor <- asv_sel_sarc_only_tocor[,apply(asv_sel_sarc_only_tocor, 2, function(x){sum(x>0)})>10]

dim(asv_sel_sarc_only)
asv_sel_sarc_only_tocor[,1:50]
# asv_sel <- asv_sel[,apply(asv_sel, 2, function(x){sum(x>0)})>50]
clin_sig_data_sarc_only_num <- clin_sig_data_sarc_only[,-c(1,3,9:11)]
plot(asv_sel_sarc_only$ASV746,clin_sig_data_sarc_only$t121)
plot(asv_sel_norm_only$ASV75,clin_sig_data_norm_only$BMI)
pheatmap(norm_r,display_numbers = norm_pstar,border_color = NA,cellwidth = 10,cellheight = 10,
         fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = F,cluster_cols = F,main = 'normal')
pheatmap(sarc_r,display_numbers = sarc_psatr,border_color = NA,cellwidth = 10,cellheight = 10,
         fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = F,cluster_cols = F,main = 'sarcopenia')
library(matrixStats)
library(stringr)
library(psych)
library(pheatmap)
clin_asv_cor_r_sel <- clin_asv_cor_r[rowMins(as.matrix(clin_asv_cor_p))<0.05,]
sarc_sel <- clin_asv_cor_r_sel
norm_sel <- clin_asv_cor_r_sel
# clin_asv_cor_p_sel <- clin_asv_cor_p[clin_rel_ASV,]
# clin_asv_cor_p_sel <- apply(clin_asv_cor_p_sel,2,p2star)

clin_rel_ASV <- rownames(clin_asv_cor_r_sel) #临床表型相关的ASV
clin_asv_cor_p_sel <- clin_asv_cor_p[clin_rel_ASV,]
clin_asv_cor_p_sel <- apply(clin_asv_cor_p_sel,2,p2star)
pheatmap(cor(clin_sig_data_sarc_only_num))
Sarc_ASV_tax_clin_rel <- Sarc_ASV_tax[clin_rel_ASV,]
View(Sarc_ASV_tax_clin_rel)
Sarc_ASV_tax_clin_rel_name <- str_c('G',Sarc_ASV_tax_clin_rel$Genus,sep = '_')
Sarc_ASV_tax_clin_rel_name[1] <- 'O_Clostridia UCG-014'
Sarc_ASV_tax_clin_rel_name[2:4] <- str_c(Sarc_ASV_tax_clin_rel_name[2:4],1:3,sep = '_')
identical(rownames(Sarc_ASV_tax_clin_rel),clin_rel_ASV)
rownames(clin_asv_cor_r_sel) <- Sarc_ASV_tax_clin_rel_name
?pheatmap
pheatmap::pheatmap(clin_asv_cor_r_sel,display_numbers = clin_asv_cor_p_sel,border_color = NA,fontsize = 10,fontsize_number = 15,number_color = 'darkblue')
Blood_sig_cor <- corr.test(asv_blood_pare,blood_sig_data_nona,method = 'spearman')
table(Blood_sig_cor$p<0.05)
table(Blood_sig_cor$p.adj<0.05)
Blood_sig_cor_r <- Blood_sig_cor$r
Blood_sig_cor_p <- Blood_sig_cor$p
Blood_asv_cor_r_sel <- Blood_sig_cor_r[rowMins(as.matrix(Blood_sig_cor_p))<0.01,]
Blood_rel_ASV <- rownames(Blood_asv_cor_r_sel) # ASV sig
Sarc_ASV_tax_blood_rel <- Sarc_ASV_tax[Blood_rel_ASV,]
Sarc_ASV_tax_blood_rel_name <- str_c('G',Sarc_ASV_tax_blood_rel$Genus,sep = '_')

clin_asv_cor_p_sel <- Blood_sig_cor_p[clin_rel_ASV,]
clin_asv_cor_p_sel <- apply(clin_asv_cor_p_sel,2,p2star)
library(ggplot2)
library(ggtree)
clin_asv_cor_r <- na.omit(clin_asv_cor_r)
rf_best12_model$terms
blood_comp <- melt(blood_sig_data_nona_T1, id  = 'group')
head(blood_comp)
p2 <- ggboxplot(blood_comp, x = 'group', y = 'value',
                color  = 'group',add = "jitter")

facet(p2, facet.by = "variable",scales="free_y",ncol = 3)+
  stat_compare_means(aes(group = group),label = "p.format",
                     vjust = 4.5, method = 'wilcox.test', size = 4)+
  scale_color_brewer(palette = "Set1")+
  theme_classic()
blood_index <- data.frame(
  blood_pid = full_blood_id,
  group = blood_sig_data_nona_T1$group
)

full_blood_id <- rownames(blood_sig_data_nona_T1)
full_asv_blood <- sel_ASV_box2[,full_blood_id]
full_asv_blood_tocor <- t(full_asv_blood)
blood_asv_cor <- corr.test(full_asv_blood_tocor,blood_sig_data_nona ,method = 'spearman',adjust = 'none')
blood_asv_cor_r <- blood_asv_cor$r
pheatmap(blood_asv_cor_r)
blood_asv_cor_p <- blood_asv_cor$    
  blood_asv_cor_p_sig <- blood_asv_cor_p[rowMin(blood_asv_cor_p)<0.05,]
blood_asv_cor_p_star <- apply(blood_asv_cor_p_sig,2,p2star)
blood_asv_cor_r_sig <- blood_asv_cor_r[rownames(blood_asv_cor_p_sig),]
blood_asv_cor_sig_id <- rownames(blood_asv_cor_r_sig)
pheatmap(blood_asv_cor_r_sig,display_numbers = blood_asv_cor_p_star,border_color = NA,cellwidth = 10,cellheight = 10,
         fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = F,cluster_cols = F,main = 'blood cor')
# blood group ====
blood_norm <- blood_sig_data_nona[blood_sig_data_nona_T1$group == 'Normal',]
blood_sarc <- blood_sig_data_nona[blood_sig_data_nona_T1$group != 'Normal',]
full_asv_blood_tocor_norm <- full_asv_blood_tocor[rownames(blood_norm),]
full_asv_blood_tocor_sarc <- full_asv_blood_tocor[rownames(blood_sarc),]
blood_norm_cor <- corr.test(full_asv_blood_tocor_norm,blood_norm ,method = 'spearman',adjust = 'none')
blood_norm_cor_r <- na.omit(blood_norm_cor$r)
blood_norm_cor_p <- na.omit(blood_norm_cor$p)
blood_norm_cor_p_sig <- blood_norm_cor_p[rowMin(blood_norm_cor_p)<0.05,]
blood_norm_cor_star <- apply(blood_norm_cor_p_sig,2,p2star)
blood_norm_cor_r_sig <- blood_norm_cor_r[rownames(blood_norm_cor_p_sig),]
pheatmap(blood_norm_cor_r_sig,display_numbers = blood_norm_cor_star,border_color = NA,cellwidth = 10,cellheight = 10,
         fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = F,cluster_cols = F,main = 'blood cor norm')

blood_sarc_cor <- corr.test(full_asv_blood_tocor_sarc,blood_sarc ,method = 'spearman',adjust = 'none')

blood_sarc_cor_r <- na.omit(blood_sarc_cor$r)
blood_sarc_cor_p <- na.omit(blood_sarc_cor$p)
blood_sarc_cor_p_sig <- blood_sarc_cor_p[rowMin(blood_sarc_cor_p)<0.05,]
blood_sarc_cor_star <- apply(blood_sarc_cor_p_sig,2,p2star)
blood_sarc_cor_r_sig <- blood_sarc_cor_r[rownames(blood_sarc_cor_p_sig),]
pheatmap(blood_sarc_cor_r_sig,display_numbers = blood_sarc_cor_star,border_color = NA,cellwidth = 10,cellheight = 10,
         fontsize = 10,fontsize_number = 15,number_color = 'darkblue',cluster_rows = F,cluster_cols = F,main = 'blood cor sarc')

