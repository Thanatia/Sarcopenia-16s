library(phyloseq)
library(DESeq2)
library(dplyr)
library(stringr)
colData <-  data.frame(row.names = rownames(sample_data(Sarc_ASV)), group_list= Sarc_diag) 
exprMtx <- as.data.frame(otu_table(Sarc_ASV))
exprMtx <- exprMtx[rowMeans(exprMtx)>1,]
feilinzhi <- apply(exprMtx, 1, function(x){sum(x>0)})
exprMtx <- exprMtx[feilinzhi>20,]
exprMtx <- exprMtx[order(rowMadDiffs(as.matrix(exprMtx)),decreasing = T),]
exprMtx[,exprMtx[1,]==0] <- 1
boxplot(as.numeric(exprMtx['ASV41',])~Sarc_diag)
group_list <- Sarc_diag
dds = DESeqDataSetFromMatrix(countData = as.matrix(exprMtx), colData = colData, design= ~ group_list)
dds = DESeq(dds)
res1 = as.data.frame(results(dds))
res1 <- na.omit(res1)
head(res1)
plot(res1$log2FoldChange,-log10(res1$padj))
logFC_cutoff <- 1
library(ggplot2)
gaibian <- ifelse(res1$padj<0.1,as.character(res1$change),'NOT')
table(gaibian)
ggplot(res1,aes(x  = log2FoldChange, y = -log(pvalue),color = gaibian, size = log(baseMean)))+geom_point()+theme_bw()
res1_de <- res1 %>% filter(padj<0.1) %>% arrange(desc(abs(log2FoldChange)))
TOP_ASV <- rownames(res1_de)[1:10]
pheatmap::pheatmap(exprMtx[TOP_ASV,],annotation_col = annot_col,scale = 'row')

sel_ASV_boxplot <- as.data.frame(log(t(exprMtx[TOP_ASV,])))%>%mutate(group = Sarc_diag)
colnames(sel_ASV_boxplot_hi)[2:3] <- c('ASV','log10_Counts')
p2<-ggboxplot(sel_ASV_boxplot_hi_top10, x = "group", y = 'log10_Counts',
              color = 'group',#palette = "npg",
              xlab=FALSE,
              add = "jitter",
              #outlier.shape = NA,
              #ylim=(c(0,30000))
)
facet(p2 , facet.by = "ASV",scales="free_y",ncol = 3)+
  stat_compare_means(aes(group = group),label = "p.format",
                     vjust = 1, method = 'wilcox.test', size = 4)+
  scale_color_brewer(palette = "Set1")+
  theme_classic()