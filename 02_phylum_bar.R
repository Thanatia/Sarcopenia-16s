# color board
set_colors <- c(
  "#DA5724",  "#508578", "#c481fd",  "#CD9BCD","#AD6F3B", "#652926", "#6dcff6",  "#C84248",  "#1e90ff", "#8569D5", 
  "#cc0953", "#D1A33D", "grey",  "pink", "gold","#8A7C64", "#599861", "navy", "#5F7FC7" , "tomato", "#673770",  
  "#008080", "#2F4F4F", "#FAEBD7", "#ff1493", "#5e738f","#808000", "#D14285", "#ffa500", "cbd588", "wheat", 
  "#d2b48c", "cyan2","black",  "#BC8F8F", "#800000","#008B8B",  "#BC8F8F"
)
# Removing singletons
RANK = "Phylum"
plot_bar(Sarc_ASV_fr, fill = "Phylum")
dev.off()
phylumGlommed = tax_glom(Sarc_ASV_fr, RANK)
# Sarc_ASV_prop_phy <- prune_taxa(taxa_sums(phylumGlommed) > 0, phylumGlommed)
Sarc_ASV_prop_phy <- transform_sample_counts(phylumGlommed, function(x) x/sum(x))
Sarc_ASV_prop_phy_stack <-  merge_samples(Sarc_ASV_prop_phy, "Sarcopenia")
Sarc_ASV_prop_phy_stack <- transform_sample_counts(Sarc_ASV_prop_phy_stack, function(x) x/sum(x))
plot_bar(Sarc_ASV_prop_phy_stack, fill = RANK)+scale_fill_manual(values=mycol_16)+theme_classic()
# Sarc_ASV_prop_phy <- tax_glom(Sarc_ASV_prop_phy, RANK)
# RANK <- 'Phylum'
# Sarc_ASV_prop_ <- tax_glom(GPfr, RANK)
OTU_Phylum <- as.data.frame(t(otu_table(Sarc_ASV_prop_phy)))
rowSums(OTU_Phylum)
OTU_Phylum <- as.data.frame(t(OTU_Phylum))
OTU_Phylum$group <- Diag
OTU_Phylum$Sample <- rownames(OTU_Phylum)
tax_table(Sarc_ASV_prop_phy)
OTU_Phylum <- OTU_Phylum[order(OTU_Phylum$ASV15,decreasing = T),]
OTU_Phylum$Sample <- factor(OTU_Phylum$Sample,levels = OTU_Phylum$Sample)
OTU_Phylum_bar <- melt(OTU_Phylum, by = 'group')
colnames(OTU_Phylum_bar)[4] <- 'Aboundance'
colnames(OTU_Phylum_bar)[3] <- 'ASV'
TAX_Phylum <- as.data.frame(tax_table(Sarc_ASV_prop_phy))
TAX_Phylum$ASV <- rownames(TAX_Phylum)
OTU_Phylum_bar2 <- full_join(OTU_Phylum_bar,TAX_Phylum)
head(OTU_Phylum_bar2)
ggplot(OTU_Phylum_bar2,aes(x = Sample, y = Aboundance, fill = Phylum ))+
  geom_bar(stat = "identity")+scale_fill_manual(values=mycol_16)+theme_classic()+
  facet_wrap(~group, scale = 'free')
# F/B 
tax_table(Sarc_ASV_prop_phy)
colSums(ASV_phy_data)
ASV_phy_data <- as.data.frame(t(otu_table(Sarc_ASV_prop_phy)))
boxplot(ASV_phy_data$ASV15/ASV_phy_data$ASV34~Diag)
Bacteroidota <- ASV_phy_data$ASV34
Firmicutes <- ASV_phy_data$ASV15
FB_data <- data.frame(
  FB_ratio = log(Firmicutes)-log(Bacteroidota+0.00001),
  Group = Diag
)

FB_p <-ggboxplot(FB_data, x = "Group", y = 'FB_ratio',
                 color = 'Group',
                 xlab=FALSE,
                 add = "jitter",
d
)

FB_p+theme_bw()+stat_compare_means(aes(group = Group),label = "p.format",
                                   vjust = 1, method = 'wilcox.test', size = 4)+
  scale_color_brewer(palette = "Set1")
