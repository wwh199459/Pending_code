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
A2 <- read.csv(file= "lungcancer.csv", header = T, row.names = 1)
mydata <- CreateSeuratObject(counts = A2, project = "lung") 
mydata <- NormalizeData(mydata)
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)
saveRDS(mydata, file="scRNAlungcancer.rds")

#使用小提琴图可视化QC指标
plot <- VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot
#nCount_RNA与nFeature_RNA的相关性
plot2 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

scRNA <- mydata

scRNA <- readRDS("scRNAlungcancer.rds")
##meta.data
proj_name <- data.frame(proj_name=row.names(scRNA@meta.data))
rownames(proj_name) <- row.names(scRNA@meta.data)

library(stringr)
###按_拆分###
b <- str_split_fixed(proj_name$proj_name, "_", 5)
rownames(b) <- row.names(scRNA@meta.data)
b<-as.matrix(b) 
d = data.frame(b)

###单细胞增加列####
scRNA <- AddMetaData(scRNA, d)
saveRDS(scRNA, file="scRNAlungcancerprime.rds")

scRNA <- readRDS("scRNAasthma.rds")
####group命名####
scRNA$group<-recode(scRNA$X3,
                     "PTR1" = "PT",
                     "PTR2" = "PT",
                     "PTR3" = "PT",
                     "PTR4" = "PT",
                     "PTR6" = "PT",
                     "N" = "NAT",
                     "N1" = "NAT",
                     "N2" = "NAT",
                     "M" = "RT")

saveRDS(scRNA, file="scRNAlungcancersecond.rds")

##
scRNA <- readRDS("scRNAlungcancersecond.rds")
DefaultAssay(scRNA) <- "RNA"

##
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)
col.num <- length(levels(as.factor(scRNA@meta.data$X3)))


dir.create('QC')
violin <-VlnPlot(scRNA, group.by = "X3",  
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
print(c("minGene=1000", "maxGene=12000", "pctMT=20"))
minGene=1000
maxGene=10000
pctMT=10

##
scRNA1 <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
table(scRNA1$group)
col.num <- length(levels(as.factor(scRNA@meta.data$X3)))
violin <-VlnPlot(scRNA, group.by = "X3",
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
scRNA <- readRDS("QC/scRNA.rds")
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
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="X3") 
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster1/pca.pdf", plot = plotc, width = 8, height = 4) 
ggsave("cluster1/pca.png", plot = plotc, width = 8, height = 4)

pc.num=1:20

######
scRNA <- FindNeighbors(scRNA, dims = pc.num) 
scRNA <- FindClusters(scRNA, resolution = 0.1)
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

scRNA <- readRDS("cell_identify/scRNA.rds")
scRNA$celltype<-recode(scRNA$seurat_clusters,
                       "0" = "canonical immune",
                       "1" = "canonical immune",
                       "2" = "canonical immune",
                       "3" = "Epithelial",
                       "4" = "Epithelial",
                       "5" = "Epithelial",
                       "6" = "canonical immune",
                       "7" = "Epithelial",
                       "8" = "Epithelial",
                       "9" = "Epithelial",
                       "10" = "Epithelial",
                       "11" = "canonical immune",
                       "12" = "canonical immune",
                       "13" = "stromal cells")

scRNA$celltype2<-recode(scRNA$seurat_clusters,
                       "0" = "T cells",
                       "1" = "myeloid cells",
                       "2" = "B cells",
                       "3" = "Epithelial",
                       "4" = "Epithelial",
                       "5" = "Epithelial",
                       "6" = "Fibroblasts",
                       "7" = "Epithelial",
                       "8" = "Epithelial",
                       "9" = "Epithelial",
                       "10" = "Epithelial",
                       "11" = "mast cells",
                       "12" = "Natural killer T cells",
                       "13" = "stromal cells")
saveRDS(scRNA, file="cell_identify/scRNA.rds")  

dir.create("canonical immune")
setwd("K:/mateng/canonical immune")
Cells.sub <- subset(scRNA@meta.data, celltype=="canonical immune")
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))
saveRDS(scRNAsub, file="scRNACI.rds")  

dir.create("Epithelial")
setwd("K:/mateng/Epithelial")
Cells.sub <- subset(scRNA@meta.data, celltype=="Epithelial")
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))
saveRDS(scRNAsub, file="scRNAEP.rds")\

dir.create("stromal cells")
setwd("K:/mateng/stromal cells")
Cells.sub <- subset(scRNA@meta.data, celltype=="stromal cells")
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))
saveRDS(scRNAsub, file="scRNAST.rds") 

scRNA <- readRDS("scRNAST.rds")
Cells.sub <- subset(x = scRNA, subset = ZC3H12A > 0)
table(Cells.sub$seurat_clusters)

Cells.sub$gene1<-recode(Cells.sub$X3,
                        "M" = "ZC3H12A h",
                        "PTR1" = "ZC3H12A h",
                        "PTR2" = "ZC3H12A h",
                        "PTR4" = "ZC3H12A h")

Cells.sub2 <- subset(x = scRNA, subset = ZC3H12A <= 0)
table(Cells.sub2$seurat_clusters)

Cells.sub2$gene1<-recode(Cells.sub2$X3,
                        "M" = "ZC3H12A l",
                        "N" = "ZC3H12A l",
                        "N2" = "ZC3H12A l",
                        "PTR1" = "ZC3H12A l",
                        "PTR2" = "ZC3H12A l")
scRNAST  <- merge(Cells.sub,Cells.sub2)

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
scRNA <- readRDS("cell_identify/scRNA.rds")
Cells.sub <- subset(x = scRNA, subset = ZC3H12A > 0)
Cells.sub2 <- subset(x = scRNA, subset = ZC3H12A <= 0)
table(Cells.sub$X1)
table(Cells.sub2$X1)
Cells.sub$gene1<-recode(Cells.sub$X1,
                        "SCLC" = "ZC3H12A h")
Cells.sub2$gene1<-recode(Cells.sub2$X1,
                        "SCLC" = "ZC3H12A l")
scRNAZC3H12A  <- merge(Cells.sub,Cells.sub2)
table(scRNAZC3H12A$orig.ident)
saveRDS(scRNAZC3H12A, file="cell_identify/scRNAZC3H12A.rds")

scRNA <- readRDS("cell_identify/scRNA.rds")
Cells.sub <- subset(x = scRNA, subset = IFIH1 > 0)
Cells.sub2 <- subset(x = scRNA, subset = IFIH1 <= 0)
table(Cells.sub$X1)
table(Cells.sub2$X1)
Cells.sub$gene2<-recode(Cells.sub$X1,
                        "SCLC" = "IFIH1 h")
Cells.sub2$gene2<-recode(Cells.sub2$X1,
                         "SCLC" = "IFIH1 l")
scRNAIFIH1  <- merge(Cells.sub,Cells.sub2)
table(scRNAIFIH1$gene2)
saveRDS(scRNAIFIH1, file="cell_identify/scRNAIFIH1.rds")
table(scRNAIFIH1$gene2)

#数据准备，分别创建CellChat对象
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
scRNA <- readRDS("cell_identify/scRNAZC3H12A.rds")
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
options(stringsAsFactors = FALSE)
table(scRNA$celltype2)
table(scRNA$gene1)
Idents(scRNA) <- 'gene1'
sco.L <- subset(scRNA, idents = c('ZC3H12A h'))########实验组
sco.P <- subset(scRNA, idents = c('ZC3H12A l'))########对照组

#创建cellchat对象
cco.L <- createCellChat(sco.L@assays$RNA@data, meta = sco.L@meta.data, group.by = "celltype2")
cco.P <- createCellChat(sco.P@assays$RNA@data, meta = sco.P@meta.data, group.by = "celltype2")
dir.create("vs")
setwd("./vs")
save(cco.L, cco.P, file = "cco.rda")

############分析对照组样本cco.L的细胞通讯网络
cellchat <- cco.P
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.P <- cellchat
saveRDS(cco.P, "cco.P.rds")

#########分析实验组样本cco.L的细胞通讯网络
cellchat <- cco.L
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.L <- cellchat
saveRDS(cco.L, "cco.L.rds")

#合并cellchat对象
cco.list <- list(pbmc=cco.P, til=cco.L)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

#3. 可视化</h3>
#3.1 所有细胞群总体观：通讯数量与强度对比</h5>
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)
p

#数量与强度差异网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# save as Diff_number_strength_net.pdf

#数量与强度差异热图</p>红色是case相对于control上调的，蓝色是下调的
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
# save as Diff_number_strength_heatmap.pdf

#细胞互作数量对比网络图</p>
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
# save as Counts_Compare_net.pdf

#########3.2 指定细胞互作数量对比网络图</h5>
par(mfrow = c(1,2))
s.cell <- c("B cells", "Epithelial", "Fibroblasts", "mast cells", "myeloid cells", "Natural killer T cells", "stromal cells", "T cells")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

####3.3 保守和特异性信号通路的识别与可视化</h5>
## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 15, height = 12)   

#####3.4 流行学习识别差异信号通路
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

saveRDS(cellchat, "cellchat.rds") 

##3.5 细胞信号模式对比</h5>
library(ComplexHeatmap)
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_all.pdf  10*6 

##输出信号模式对比</p>
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_outgoing.pdf  10*6 

#输入信号模式对比</p>
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_incoming.pdf  10*6

#特定信号通路的对比</h5>
#normal存在，疾病不存在
#CD40
#CD86
#LIGHT
#NECTIN
#TIGIT
#SEMA4
#网络图</p>
pathways.show <- c("MHC−I") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_net.pdf  10*6.5

#热图</p>
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# save as Compare_IL16_heatmap.pdf  12*6.5

#和弦图</p>
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_chord.pdf  10*6.5

#######配体-受体对比分析 
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7,8), targets.use = c(2),  comparison = c(1, 2), angle.x = 45)
ggsave("Compare_LR_bubble.pdf", p, width = 25, height = 30)

#气泡图展示上调或下调的配体受体对</p>
p1 <- netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7,8), targets.use = c(2), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7,8), targets.use = c(2), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Compare_LR_regulated.pdf", p1, width = 25, height = 30)

#和弦图</p>
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use = c(1,3,4,5,6,7,8), targets.use = c(2), signaling = "FN1", 
                       lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("Signaling from Treg - ", names(cco.list)[i]))
}
# save as Compare_LR_chord.pdf  10*6.5


#数据准备，分别创建CellChat对象IFIH1
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
scRNA <- readRDS("cell_identify/scRNAIFIH1.rds")
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
options(stringsAsFactors = FALSE)
table(scRNA$celltype2)
table(scRNA$gene2)
Idents(scRNA) <- 'gene2'
sco.L <- subset(scRNA, idents = c('IFIH1 h'))########实验组
sco.P <- subset(scRNA, idents = c('IFIH1 l'))########对照组

#创建cellchat对象
cco.L <- createCellChat(sco.L@assays$RNA@data, meta = sco.L@meta.data, group.by = "celltype2")
cco.P <- createCellChat(sco.P@assays$RNA@data, meta = sco.P@meta.data, group.by = "celltype2")
dir.create("vs2")
setwd("./vs2")
save(cco.L, cco.P, file = "cco.rda")

############分析对照组样本cco.L的细胞通讯网络
cellchat <- cco.P
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.P <- cellchat
saveRDS(cco.P, "cco.P.rds")

#########分析实验组样本cco.L的细胞通讯网络
cellchat <- cco.L
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.L <- cellchat
saveRDS(cco.L, "cco.L.rds")

#合并cellchat对象
cco.list <- list(pbmc=cco.P, til=cco.L)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

#3. 可视化</h3>
#3.1 所有细胞群总体观：通讯数量与强度对比</h5>
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)
p

#数量与强度差异网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# save as Diff_number_strength_net.pdf

#数量与强度差异热图</p>红色是case相对于control上调的，蓝色是下调的
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
# save as Diff_number_strength_heatmap.pdf

#细胞互作数量对比网络图</p>
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
# save as Counts_Compare_net.pdf

#########3.2 指定细胞互作数量对比网络图</h5>
par(mfrow = c(1,2))
s.cell <- c("B cells", "Epithelial", "Fibroblasts", "mast cells", "myeloid cells", "Natural killer T cells", "stromal cells", "T cells")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

####3.3 保守和特异性信号通路的识别与可视化</h5>
## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 15, height = 12)   

#####3.4 流行学习识别差异信号通路
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

saveRDS(cellchat, "cellchat.rds") 

##3.5 细胞信号模式对比</h5>
library(ComplexHeatmap)
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_all.pdf  10*6 

##输出信号模式对比</p>
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_outgoing.pdf  10*6 

#输入信号模式对比</p>
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_incoming.pdf  10*6

#特定信号通路的对比</h5>
#normal存在，疾病不存在
#CD40
#CD86
#LIGHT
#NECTIN
#TIGIT
#SEMA4
#网络图</p>
pathways.show <- c("MHC−I") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_net.pdf  10*6.5

#热图</p>
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# save as Compare_IL16_heatmap.pdf  12*6.5

#和弦图</p>
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_chord.pdf  10*6.5

#######配体-受体对比分析 
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7,8), targets.use = c(2),  comparison = c(1, 2), angle.x = 45)
ggsave("Compare_LR_bubble.pdf", p, width = 25, height = 30)

#气泡图展示上调或下调的配体受体对</p>
p1 <- netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7,8), targets.use = c(2), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7,8), targets.use = c(2), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Compare_LR_regulated.pdf", p1, width = 25, height = 30)

#和弦图</p>
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use = c(1,3,4,5,6,7,8), targets.use = c(2), signaling = "WNT", 
                       lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("Signaling from Treg - ", names(cco.list)[i]))
}
# save as Compare_LR_chord.pdf  10*6.5

















