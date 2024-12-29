rm(list = ls())
library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(patchwork)

setwd("K:/mateng")
A2 <- read.csv(file= "lungcancer.csv", header = T)
mydata <- CreateSeuratObject(counts = A2, project = "PBMC") 
mydata <- NormalizeData(mydata)
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)
saveRDS(mydata, file="scRNAasthma.rds")

#使用小提琴图可视化QC指标
plot <- VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot
#nCount_RNA与nFeature_RNA的相关性
plot2 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

scRNA <- mydata

scRNA <- readRDS("scRNAasthma.rds")
##meta.data
proj_name <- data.frame(proj_name=row.names(scRNA@meta.data))
rownames(proj_name) <- row.names(scRNA@meta.data)

library(stringr)
###按_拆分###
b <- str_split_fixed(proj_name$proj_name, "_", 2)
rownames(b) <- row.names(scRNA@meta.data)
b<-as.matrix(b) 
d = data.frame(b)

###单细胞增加列####
scRNA <- AddMetaData(scRNA, d)
saveRDS(scRNA, file="scRNAasthma.rds")

scRNA <- readRDS("scRNAasthma.rds")
####group命名####
scRNA$sample<-recode(scRNA$X2,
                       "350C" = "SA1",
                       "865C" = "SA2",
                       "866C" = "SA3",
                       "A307C" = "SA4",
                       "A311C" = "SA5",
                       "350P" = "HC1",
                       "865P" = "HC2",
                       "866P" = "HC3",
                       "A307P" = "SAPOL1",
                       "A311P" = "SAPOL2",
                       "459C" = "SAPOL3",
                       "851C" = "SAPOL4",
                       "868C" = "SAPOL5",
                       "459P" = "HCPOL1",
                       "851P" = "HCPOL2",
                       "868P" = "HCPOL3")

saveRDS(scRNA, file="scRNAasthma.rds")

##
DefaultAssay(scRNA) <- "RNA"

##
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)
col.num <- length(levels(as.factor(scRNA@meta.data$sample)))


dir.create('QC')
violin <-VlnPlot(scRNA, group.by = "sample",  
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.01, #
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)

##QC###
print(c("minGene=200", "maxGene=6000", "pctMT=10"))
minGene=200
maxGene=6000
pctMT=10

##
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(as.factor(scRNA@meta.data$sample)))
violin <-VlnPlot(scRNA, group.by = "sample",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 6)

##
saveRDS(scRNA, file="QC/scRNA.rds")


########寻找高变基因###
dir.create("cluster")
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
ggsave("cluster/VariableFeatures.pdf", plot = plot, width = 8, height = 6) 
ggsave("cluster/VariableFeatures.png", plot = plot, width = 8, height = 6)

##
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)

#
GetAssayData(scRNA,slot="counts",assay="RNA") 
#             
GetAssayData(scRNA,slot="data",assay="RNA")
#
GetAssayData(scRNA,slot="scale.data",assay="RNA") 
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA))
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)

######
dir.create('cluster1')
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="sample") 
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster1/pca.pdf", plot = plotc, width = 8, height = 4) 
ggsave("cluster1/pca.png", plot = plotc, width = 8, height = 4)

pc.num=1:20

######
scRNA <- FindNeighbors(scRNA, dims = pc.num) 
scRNA <- FindClusters(scRNA, resolution = 0.4)
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster1/cell_cluster.csv',row.names = F)

##########
#tSNE
scRNA = RunTSNE(scRNA, dims = pc.num)
embed_tsne <- Embeddings(scRNA, 'tsne')
write.csv(embed_tsne,'cluster1/embed_tsne.csv')
plot1 = DimPlot(scRNA, reduction = "tsne") 
ggsave("cluster1/tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("cluster1/tSNE.png", plot = plot1, width = 8, height = 7)
#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap,'cluster1/embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap") 
ggsave("cluster1/UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("cluster1/UMAP.png", plot = plot2, width = 8, height = 7)
#
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("cluster1/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
ggsave("cluster1/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
##
saveRDS(scRNA, file="cluster1/scRNAafter.rds")

########Cluster
dir.create("cell_identify")
#
diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "cell_identify/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "cell_identify/top10_diff_genes_wilcox.csv", row.names = F)
#
diff.mast = FindAllMarkers(scRNA, test.use = 'MAST')
all.markers = diff.mast %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "cell_identify/diff_genes_mast.csv", row.names = F)
write.csv(top10, "cell_identify/top10_diff_genes_mast.csv", row.names = F)
#bulkRNA
diff.deseq2 = FindAllMarkers(scRNA, test.use = 'DESeq2', slot = 'counts')
all.markers = diff.deseq2 %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "cell_identify/diff_genes_deseq2.csv", row.names = F)
write.csv(top10, "cell_identify/top10_diff_genes_deseq2.csv", row.names = F)

##top10
top10_genes <- read.csv("cell_identify/top10_diff_genes_wilcox.csv")
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNA)) 
plot1 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave("cell_identify/top10_markers.pdf", plot=plot1, width=8, height=6) 
ggsave("cell_identify/top10_markers.png", plot=plot1, width=8, height=6)
saveRDS(scRNA, file="cell_identify/scRNA.rds")

scRNA <- readRDS("cell_identify/scRNA.rds")
scRNA$celltype<-recode(scRNA$seurat_clusters,
                       "0" = "CD4+ T cells",
                       "1" = "CD8+ T cells",
                       "2" = "NK cells",
                       "3" = "CD4+ T cells",
                       "4" = "B cells",
                       "5" = "Monocytes")
saveRDS(scRNA, file="cell_identify/scRNA.rds")  

dir.create("cell_identify5top")
#挑选部分基因
select_genes <- c('LDHB',	'TCF7',	'RCAN3',	'LEF1',	'RPS3AP6',	'PASK',	'LRRN3',	'RP4-594I10.3',	'NELL2',	'LDLRAP1')
#vlnplot展示
p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="celltype", ncol=5)
ggsave("cell_identify5top/selectgenes_VlnPlot_0.png", p1, width=16,height=8)

#挑选部分基因
select_genes <- c('KLRD1',	'NKG7',	'GNLY',	'FGFBP2',	'GZMH',	'CST7',	'GZMA',	'GZMB',	'C1orf21',	'PRF1')
#vlnplot展示
p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="celltype", ncol=5)
ggsave("cell_identify5top/selectgenes_VlnPlot_1.png", p1, width=16,height=8)

#挑选部分基因
select_genes <- c('IFI44L',	'IFI6',	'IFI44',	'IFIT3',	'MX1',	'IFIT1',	'MX2',	'PARP9',	'RSAD2',	'HERC6')
#vlnplot展示
p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="celltype", ncol=5)
ggsave("cell_identify5top/selectgenes_VlnPlot_2.png", p1, width=16,height=8)

#挑选部分基因
select_genes <- c('IGHM',	'MS4A1',	'CD79A',	'BANK1',	'IGHD',	'HLA-DRA',	'MEF2C',	'IGLC2',	'HLA-DMB',	'LINC00926')
#vlnplot展示
p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="celltype", ncol=5)
ggsave("cell_identify5top/selectgenes_VlnPlot_4.png", p1, width=16,height=8)

#挑选部分基因
select_genes <- c('LYZ',	'VMO1',	'MS4A7',	'AIF1',	'CYBB',	'IGSF6',	'LST1',	'NCF2',	'FGL2',	'CST3')
#vlnplot展示
p1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="celltype", ncol=5)
ggsave("cell_identify5top/selectgenes_VlnPlot_5.png", p1, width=16,height=8)


library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
rm(list=ls())
setwd("D:/asthmawen/GSE172495_RNA_Matrix_PBMC.csv")
dir.create("CD4")
scRNA <- readRDS("cell_identify/scRNA.rds")
##提取细胞子集
Cells.sub <- subset(scRNA@meta.data, celltype=="T cells")
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))






















