#REQUIRED PACKAGES 
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
remotes::install_github(repo = 'mojaveazure/seurat-disk', ref = 'develop')
remotes::install_github("FASTGenomics/fgread-r")
install.packages("Seurat")
devtools::install_github("immunogenomics/harmony") #if you have problem installing it, see: https://github.com/immunogenomics/symphony/issues/2
BiocManager::install("celldex")

#LIBRARIES
library(fgread)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(celldex)
library(dplyr)
library(scales)


#Import source data (SEURAT VISUALIZATION) for heatmap and dotplot clustering and ordering  
source("/Users/otavio.cmarques/Documents/Otavio/Orientacoes /Lena/Covid Datasets/10x PBMC/Usage/Scripts/SEURAT VISUALISATION.R")

#Import Data Set
# load data from rds file: folder address + seurat file name
candida <- readRDS("/Users/otavio.cmarques/Documents/Otavio/Orientacoes /Ranieri/scRNASEQ/seurat_object_final_anonymized.Rds")

# uptade seurat object  
candida <- UpdateSeuratObject(object = candida)

View(candida@meta.data)

#View(candida@assays$RNA@scale.data)
levels(candida)

## Normalization
candida <- NormalizeData(object = candida, normalization.method = "LogNormalize", scale.factor = 1e4)
## Define variable genes
candida <- FindVariableFeatures(object = candida, assay="RNA", selection.method = 'vst', nfeatures = 10000)
## Scaling
candida <- ScaleData(object = candida, vars.to.regress = c("nCount_RNA"))
candida <- ScaleData(object = candida, features = rownames(candida), vars.to.regress = c("nCount_RNA"))

## PCA
candida <- RunPCA(object = candida)

## UMAP
candida <- RunUMAP(candida, reduction = "pca", dims = 1:20, seed.use = 42)

## Harmony
candida <- RunHarmony(candida, group.by.vars = "orig.ident", reduction = "pca", dims.use = 1:20)

## Harmony UMAP
candida <- RunUMAP(candida, reduction = "harmony", dims = 1:20, seed.use = 42)

# Fing neighbours
#candida <- FindNeighbors(object = candida, dims = 1:20, reduction="harmony", force.recalc = TRUE)

# Calculate clusters
#candida <- FindClusters(object = candida, resolution = 0.4, algorithm = 1)

#Rename Idents
new.cluster.ids <- c("B", "CD4T", "CD8T", "Monocyte", "NK", "pDC")
names(new.cluster.ids) <- levels(candida)
candida <- RenameIdents(candida, new.cluster.ids)

candida@meta.data$stim[candida@meta.data$stim=="Unstim"] <- "Resting"
candida@meta.data$stim[candida@meta.data$stim=="Stim"] <- "C.albicans"

#Dimplot by cells
DimPlot(object = candida, 
        reduction = 'umap',
        #group.by="active.ident",
        label = TRUE, repel = FALSE) + ggtitle("Cell clusters")+
theme(axis.text.x = element_text(size = 17),
      axis.text.y = element_text(size = 17),  
      axis.title.x = element_text(size = 17),
      axis.title.y = element_text(size = 17),
      title =element_text(size=17, face='bold'),
      plot.title = element_text(hjust = 0.5))+NoLegend()

#Dimplot divided by group
colorResting <- c("grey", "#1657E1")
colorC.albicans <- c("grey", "#ffa000")

DimPlot(candida,group.by="stim",cols=colorResting,order=c("Resting"))+NoLegend()+ggtitle("Resting")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))
DimPlot(candida,group.by="stim",cols=colorC.albicans,order=c("C.albicans"))+NoLegend()+ggtitle("C. albicans")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'),
        plot.title = element_text(face='bold.italic'))

#merge
DimPlot(object = candida, 
        reduction = 'umap',
        label = F,
        group.by = "stim",
        cols=c("#ffa000","#1657E1"))+ggtitle("Merge")+
  scale_color_manual("Group", values = c(C.albicans = "#ffa000", Resting= "#1657E1"))+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))+NoLegend()

#Cluster Marker Genes: Finding DEGs
markers <- FindAllMarkers(object = candida, test.use = "wilcox")
#topgenes <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#export your table of DEG
write.csv(markers,"DEGs.csv")

#Dotplot
features=c("ENSG00000121060",	"ENSG00000184979",	"ENSG00000134321",	"ENSG00000055332",	"ENSG00000157601",	"ENSG00000216490",	"ENSG00000140968",	"ENSG00000135114",	"ENSG00000111537",	"ENSG00000089127",	"ENSG00000090339",	"ENSG00000187608",	"ENSG00000111331",	"ENSG00000119917",	"ENSG00000119922",	"ENSG00000185745",	"ENSG00000174125",	"ENSG00000117676",	"ENSG00000170458",	"ENSG00000108861",	"ENSG00000104312",	"ENSG00000170345",	"ENSG00000107968",	"ENSG00000100906",	"ENSG00000134070",	"ENSG00000109320",	"ENSG00000111335",	"ENSG00000132274",	"ENSG00000132530",	"ENSG00000172183",	"ENSG00000162654",	"ENSG00000125347",	"ENSG00000154451",	"ENSG00000034152")
levels(candida) <- c("CD4T", "CD8T", "NK",  "B", "Monocyte","pDC")

DotPlot(candida,features = features, cols=c("#ffa000", "#1657E1"), split.by = "stim")+
  #scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n =11, name = "GnBu")))+
  theme(axis.title = element_blank())+ coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#Heatmap
cols.use <- list(stim=c(Resting= "#1657e1", C.albicans="#eaac43"))

levels(candida)
levels(candida) <- c("CD4T", "CD8T", "NK",  "B", "Monocyte","pDC")

DoMultiBarHeatmap(candida, features=features, additional.group.by = "stim", additional.group.sort.by = "stim", size = 0, cols.use = cols.use, draw.lines = TRUE, lines.width = NULL)+
  theme(text = element_text(size = 17, color = "black"))+
  scale_fill_gradient2(low="olivedrab1", mid="grey", high="royalblue4", na.value = "white")
##### DIVIDED BY GROUPS FOR THE GENE UMAP#####

#RESTING
#subset
candida_Resting <- subset(candida, subset=stim %in% c("Resting"))
candida_Resting <- subset(candida_Resting)

##C.ALBICANS
#subset
candida_C.albicans <- subset(candida, subset=stim %in% c("C.albicans"))
candida_C.albicans <- subset(candida_C.albicans)

#UMAPs of shared genes by TLRs and IFN signaling pathways

#ISG15
FeaturePlot(candida_Resting, features=c("ENSG00000187608"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1', breaks=c(0,6),limits=c(0,6), oob=squish)+ ggtitle("ISG15")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

FeaturePlot(candida_C.albicans, features=c("ENSG00000187608"), 
            order = T, ncol = 1, cols = c('lightgrey','#ffa000'))+ 
  scale_color_gradient(low='lightgrey',high = '#ffa000', breaks=c(0,6),limits=c(0,6), oob=squish)+ ggtitle("ISG15")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))


#IRF7
FeaturePlot(candida_Resting, features=c("ENSG00000185507"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1', breaks=c(0,3),limits=c(0,3), oob=squish)+ ggtitle("IRF7")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

FeaturePlot(candida_C.albicans, features=c("ENSG00000185507"), 
            order = T, ncol = 1, cols = c('lightgrey','#ffa000'))+ 
  scale_color_gradient(low='lightgrey',high = '#ffa000', breaks=c(0,3),limits=c(0,3), oob=squish)+ ggtitle("IRF7")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))


#UBC
FeaturePlot(candida_Resting, features=c("ENSG00000150991"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1', breaks=c(0,4),limits=c(0,4), oob=squish)+ ggtitle("UBC")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

FeaturePlot(candida_C.albicans, features=c("ENSG00000150991"), 
            order = T, ncol = 1, cols = c('lightgrey','#ffa000'))+ 
  scale_color_gradient(low='lightgrey',high = '#ffa000', breaks=c(0,4),limits=c(0,4), oob=squish)+ ggtitle("UBC")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))


#UBE2N
FeaturePlot(candida_Resting, features=c("ENSG00000177889"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1', breaks=c(0,3),limits=c(0,3), oob=squish)+ ggtitle("UBE2N")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

FeaturePlot(candida_C.albicans, features=c("ENSG00000177889"), 
            order = T, ncol = 1, cols = c('lightgrey','#ffa000'))+ 
  scale_color_gradient(low='lightgrey',high = '#ffa000', breaks=c(0,3),limits=c(0,3), oob=squish)+ ggtitle("UBE2N")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

#RPS27A
FeaturePlot(candida_Resting, features=c("ENSG00000143947"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1', breaks=c(0,5),limits=c(0,5), oob=squish)+ ggtitle("RPS27A")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

FeaturePlot(candida_C.albicans, features=c("ENSG00000143947"), 
            order = T, ncol = 1, cols = c('lightgrey','#ffa000'))+ 
  scale_color_gradient(low='lightgrey',high = '#ffa000', breaks=c(0,5),limits=c(0,5), oob=squish)+ ggtitle("RPS27A")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

#UBA52
FeaturePlot(candida_Resting, features=c("ENSG00000221983"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1', breaks=c(0,4),limits=c(0,4), oob=squish)+ ggtitle("UBA52")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

FeaturePlot(candida_C.albicans, features=c("ENSG00000221983"), 
            order = T, ncol = 1, cols = c('lightgrey','#ffa000'))+ 
  scale_color_gradient(low='lightgrey',high = '#ffa000', breaks=c(0,4),limits=c(0,4), oob=squish)+ ggtitle("UBA52")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))


#UBB
FeaturePlot(candida_Resting, features=c("ENSG00000170315"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1', breaks=c(0,4),limits=c(0,4), oob=squish)+ ggtitle("UBB")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

FeaturePlot(candida_C.albicans, features=c("ENSG00000170315"), 
            order = T, ncol = 1, cols = c('lightgrey','#ffa000'))+ 
  scale_color_gradient(low='lightgrey',high = '#ffa000', breaks=c(0,4),limits=c(0,4), oob=squish)+ ggtitle("UBB")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

#IRF3
FeaturePlot(candida_Resting, features=c("ENSG00000126456"), 
            order = T, ncol = 1, cols = c('lightgrey','#1657e1'))+ 
  scale_color_gradient(low='lightgrey',high = '#1657e1', breaks=c(0,3),limits=c(0,3), oob=squish)+ ggtitle("IRF3")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

FeaturePlot(candida_C.albicans, features=c("ENSG00000126456"), 
            order = T, ncol = 1, cols = c('lightgrey','#ffa000'))+ 
  scale_color_gradient(low='lightgrey',high = '#ffa000', breaks=c(0,3),limits=c(0,3), oob=squish)+ ggtitle("IRF3")+
  theme(axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),  
        axis.title.x = element_text(size = 17),
        axis.title.y = element_text(size = 17),
        title =element_text(size=17, face='bold'))

#BOXPLOTS

library(plyr)
library(dplyr)

a <- candida@meta.data
tabela <- candida@assays$RNA@scale.data
genes <- c("ENSG00000150991",	"ENSG00000187608", "ENSG00000177889",	"ENSG00000143947",	"ENSG00000221983",	"ENSG00000170315","ENSG00000185507",	"ENSG00000126456")
tabela <- tabela[genes,]
tabela <- t(tabela)
tabela<-as.data.frame(tabela)

UBC <- tabela$ENSG00000150991
ISG15 <- tabela$ENSG00000187608
IRF7 <- tabela$ENSG00000185507
UBA52 <- tabela$ENSG00000221983
UBB <- tabela$ENSG00000170315
IRF3 <- tabela$ENSG00000126456

genesID <- c("ISG15",	"IRF7",	"IRF3",	"UBB","UBC",	"UBA52")

a <- cbind(a, UBC,	ISG15,	IRF7,	UBA52,	UBB, IRF3)
df <- aggregate(a[,genesID], list(a$sample, a$stim), mean)
df <- dplyr::rename(df, sample = Group.1)
df <- dplyr::rename(df, stim = Group.2)

library(reshape2)
dfmelted <- melt(df)

my_comparisons <- list( c("Resting", "C.albicans") )
library(lemon)
library(ggpubr)

dfmelted$variable <- factor(dfmelted$variable, levels = genesID)

ggplot(dfmelted, aes(x = stim, y = value, fill= stim))+
  facet_rep_wrap(.~variable, ncol=3)+
  geom_boxplot(outlier.shape = NA, lwd = 0.2)+
  ylim(-1, 1.5)+
  stat_compare_means(method= "t.test", label = "p.signif", hide.ns = T, label.y = 1.3,
                     aes(group = stim), comparisons = my_comparisons, vjust = 0.5)+
  geom_jitter(size=0.6, position = position_jitter(.1), alpha=0.5)+
  scale_fill_manual("Group", values = c(Resting ="#1657E1", C.albicans="#ffa000"))+
  labs(x = "Genes", y= "Average Expression")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line.x=element_line(),
        strip.background = element_rect(color="white", fill="white"),
        strip.placement = "outside",
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.5,"line"),
        legend.background = element_blank(),
        axis.line.y = element_line(),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.text  = element_text(size = 11, color = "black"))+
  scale_y_continuous(breaks=seq(-1, 1.5, 0.5), limits=c(-1, 1.5))+NoLegend()






