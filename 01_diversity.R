# alpha diversity ====
Sarc_OTU1 <- t(otu_table(Sarc_ASV))
Sarc_meta <- sample_data(Sarc_ASV)
OB.OTU <- apply(Sarc_OTU1,1,function(x){sum(x>0)})
dim(Sarc_OTU1)
library(vegan)
data_richness <- estimateR(Sarc_OTU1)# calculate richness and Chao1 using vegan package
data_evenness <- vegan::diversity(Sarc_OTU1) / log(specnumber(Sarc_OTU1))               # calculate evenness index using vegan package
data_shannon <- vegan::diversity(Sarc_OTU1, index = "shannon")                        # calculate Shannon index using vegan package
data_alphadiv <- cbind( Sarc_meta$Diag, t(data_richness), data_shannon, data_evenness) # combine all indices in one data table
library(tidyverse)
library(tidyr)
library(reshape2)
head(data_alphadiv)
colnames(data_alphadiv)[1] <- 'Diagnosis'
data_alphadiv <- as.data.frame(apply(data_alphadiv, 2, as.numeric))
rownames(data_alphadiv) <- rownames(Sarc_meta)
data_alphadiv$Diagnosis <- ifelse(data_alphadiv$Diag==1,'Sarcopenia','Control')
head(data_alphadiv)
plot(data_alphadiv$S.chao1,data_alphadiv$S.ACE)
data_alphadiv_sig <- data_alphadiv[,-c(4,7,8)]
data_alphadiv_nosig <- data_alphadiv[,c(1,4,7,8)]
# data_alphadiv$observed.otu <- OB.OTU
dd <- melt(data_alphadiv_sig,id = 'Diagnosis')
dd2 <-  melt(data_alphadiv_nosig,id = 'Diagnosis')
# dd <- dd[-(1:236),]
p2_2<-ggboxplot(dd2, x = "Diagnosis", y = 'value',
                color = 'Diagnosis',
                xlab=FALSE,
                add = "jitter",
)

facet(p2_2 , facet.by = "variable",scales="free_y")+
  stat_compare_means(aes(group = Diagnosis),label = "p.format",
                     vjust = 1, method = 't.test', size = 4)+
  scale_color_brewer(palette = "Set1")+
  theme_classic()
# beta diversity =====
# Beta diversity===
library(factoextra)
library(FactoMineR)
library(phyloseq)
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
for (i in 1:7) {
  plot_ordination(Sarc_ASV, ordinate(Sarc_ASV,method= ord_meths[i], distance=dist),
                  "samples", color="Diag")+theme_bw()+
    scale_color_manual(values = c('blue','red'))+
    ggtitle(ord_meths[i])
}

p2+stat_compare_means()
# PLS-DA ====
library(mixOmics)
Sarc_diag <- sample_data(Sarc_ASV_fr)$Sarcopenia
plsda.data <- plsda(t(otu_table(Sarc_ASV)), Sarc_diag,ncomp = 2)
plotIndiv(plsda.data,ind.names = FALSE,ellipse = TRUE,
          legend = TRUE,legend.title = "Group",
          title = "PLS-DA",X.label = "PLS1",Y.label = "PLS2",
          cex = 1,size.xlabel = 14,size.ylabel = 14,
          size.axis = 14,size.title = 14,pch = 16,
          col.per.group = c("#56B4E9", "#ff5a5a"),
          centroid = T,guide = "none")+theme_classic()
plsda_x <- plsda.data$variates$X
plot(plsda_x)
min(plsda_x)
# beta anosim=====
library(vegan)
anosim_val = anosim(plsda_x, factor(Sarc_diag), permutations=999)
anosim_val
plot(anosim_val, col = c('#5ece46',"#56B4E9", "#ff5a5a"), xlab = "Diagnosis",ylab = "Distance")