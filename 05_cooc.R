library(igraph)
library(psych)
ASV_table_data <- as.data.frame(asv_table(Sarc_ASV ))
dim(asv_table_data)
Normal_asv <- ASV_table_data[,Sarc_diag == 'Normal']
Sarc_asv <- ASV_table_data[,Sarc_diag == 'Sarcopenia']
feiling <- function(x){sum(x!=0)}
norm_feiling <- apply(Normal_asv, 1, feiling)
Normal_asv_hi <- Normal_asv[norm_feiling>10,]
sarc_feiling <- apply(Sarc_asv, 1, feiling)
Sarc_asv_hi <- Sarc_asv[sarc_feiling>10,]
Norm_cor = corr.test(t(Normal_asv_hi),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
Norm_cor_r <- Norm_cor$r
Norm_cor_p <- Norm_cor$p.adj
table(abs(Norm_cor_r)>0.5)
Norm_sel_r <- Norm_cor_r
Norm_sel_r[Norm_cor_p > 0.05] <- 0
Norm_sel_r[abs(Norm_sel_r)<0.6] <- 0
table(Norm_sel_r!=0)
Norm_igraph = graph_from_adjacency_matrix(Norm_sel_r,mode="undirected",weighted=TRUE,diag=FALSE)
Norm_igraph
Sarc_cor <-  corr.test(t(Sarc_asv_hi),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
Sarc_cor_r <- Sarc_cor$r
Sarc_cor_p <- Sarc_cor$p.adj
table(abs(Sarc_cor_r)>0.5)
Sarc_sel_r <- Sarc_cor_r
Sarc_sel_r[Sarc_cor_p > 0.05] <- 0
Sarc_sel_r[abs(Sarc_sel_r)<0.6] <- 0
table(Sarc_sel_r!=0)
Sarc_igraph = graph_from_adjacency_matrix(Sarc_sel_r,mode="undirected",weighted=TRUE,diag=FALSE)
Sarc_igraph
#Normal igraph====
bad.vsn = V(Norm_igraph)[igraph::degree(Norm_igraph) == 0]
Norm_igraph = delete.vertices(Norm_igraph, bad.vsn)
Norm_igraph
Norm_igraph.weight = E(Norm_igraph)$weight
E(Norm_igraph)$weight = NA
set.seed(123)
plot(Norm_igraph,main="Normal Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
E.color = Norm_igraph.weight
E.color = ifelse(E.color>0, "red",ifelse(E.color<0, "blue","grey"))
E(Norm_igraph)$color = as.character(E.color)
set.seed(123)
plot(Norm_igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
E(Norm_igraph)$width = abs(Norm_igraph.weight)*4
set.seed(123)
plot(Norm_igraph,main="Normal Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

# Sarcopenia igraph===
Sarc_igraph

bad.vss = V(Sarc_igraph)[igraph::degree(Sarc_igraph) == 0]
Sarc_igraph = delete.vertices(Sarc_igraph, bad.vss)
Sarc_igraph

Sarc_igraph.weight = E(Sarc_igraph)$weight

E(Sarc_igraph)$weight = NA

set.seed(123)
plot(Sarc_igraph,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

sum(Sarc_igraph.weight>0)# number of postive correlation
sum(Sarc_igraph.weight<0)# number of negative correlation

E.color = Sarc_igraph.weight
E.color = ifelse(E.color>0, "red",ifelse(E.color<0, "blue","grey"))
E(Sarc_igraph)$color = as.character(E.color)

set.seed(123)
plot(Sarc_igraph,main="Sarcopenia Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

E(Sarc_igraph)$width = abs(Sarc_igraph.weight)*4

set.seed(123)
plot(Sarc_igraph,main="Sarcopenia Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,
     vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))

Sarc_tax <- as.data.frame(tax_table(Sarc_ASV))
igraph.col = Sarc_tax[V(igraph)$name,]
unique(igraph.col$Phylum)
igraph.col$Phylum2 <- factor(igraph.col$Phylum)
table(igraph.col$Phylum2)
phylum_color <- c( 'pink', "lightblue","purple","#FF7892",'orange','#46cc4f','green','brown','grey60')
names(phylum_color) <- unique(igraph.col$Phylum2)
levels(igraph.col$Phylum2) = phylum_color #
V(igraph)$color = as.character(igraph.col$Phylum2)
color_anno <- table(V(igraph)$color,igraph.col$Phylum2)
set.seed(123)
V(igraph)$color

plot(igraph,main="Norm Co-occurrence network",vertex.frame.color=NA,
     vertex.label = NA,
     edge.lty = 1 ,edge.curved=TRUE,margin=c(0,0,0,0))

igraph_normal_for_export <- igraph
igraph <- Sarc_igraph
igraph.size = asv_table_data[V(igraph)$name,] #
igraph.size1 = rowMeans(igraph.size) 
V(igraph)$size = log2(igraph.size1+1)+3
Sarc_tax <- as.data.frame(tax_table(Sarc_ASV))
igraph.col = Sarc_tax[V(igraph)$name,]
rownames(igraph.col)[is.na(igraph.col$Phylum)]
unique(igraph.col$Phylum)
igraph.col$Phylum <- factor(igraph.col$Phylum)
table(igraph.col$Phylum)
color_anno
cc1 <- igraph.col$Phylum
p1_anno <- data.frame(
  phylum = cc1,
  phylum_color = igraph.col$Phylum
)
p1_anno <- p1_anno[!duplicated(p1_anno$phylum),]
fix(p1_anno)
save(p1_anno,file = 'color_anno.Rdata')
phylum_color <- c( 'pink', "lightblue","yellow","#FF7892",'#46cc4f','green','brown','grey60')
names(phylum_color) <- unique(igraph.col$Phylum)
levels(igraph.col$Phylum) = phylum_color 
V(igraph)$color = as.character(igraph.col$Phylum)
table(V(igraph)$color,igraph.col$Phylum)
set.seed(123)
V(igraph)$color
plot(igraph,main="Sarc Co-occurrence network",vertex.frame.color=NA,
     vertex.label=NA,
     edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))
igraph_sarc_for_export <- igraph
legend("bottomright",legend= unique(igraph.col$Phylum),col = as.character(phylum_color) ,lwd = 2) # 画出legend框

# BiocManager::install("RCy3") export to cytoscape
library(RCy3)
cytoscapePing()
createNetworkFromIgraph(Norm_igraph,"Norm_Igraph_goodgood")
createNetworkFromIgraph(igraph_sarc_for_export,"Sarc_Igraph_goodgood")
library(venn)
Sarcopenia_ASV <- unique(V(Sarc_igraph)$name)
Normal_ASV <- unique(V(Norm_igraph)$name)
Sarcopenia_only <- setdiff(Sarcopenia_ASV,Normal_ASV)
Shared_ASV <- intersect(Sarcopenia_ASV,Normal_ASV)
Normal_only <-  setdiff(Normal_ASV,Sarcopenia_ASV)
a1<-list(Normal_ASV,Sarcopenia_ASV)
names(a1) <- c('Normal','Sarcopenia')
venn(a1,zcolor=c('blue','red'),ggplot = T)