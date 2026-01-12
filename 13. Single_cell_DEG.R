library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
library(harmony)
library(readxl)



dat1 <- readRDS('single_cell_dat.rds')
obj1 <- subset(dat1,subclass == 'Immune')


obj2 <- CreateSeuratObject(counts = obj1@assays$RNA@counts,
                           meta.data = obj1@meta.data)
obj2 <- NormalizeData(obj2)
obj2 <- FindVariableFeatures(obj2)
obj2 <- ScaleData(obj2)
obj2 <- RunPCA(obj2)
obj2 <- RunUMAP(obj2,dims = 1:15)




CD8_obj <- subset(obj2,seurat_clusters %in% c(2,4,5))
CD8_obj$group1 <- 0
CD8_obj$group1[CD8_obj$group=='DKD'] <- 1
CD8_marker <- FindMarkers(CD8_obj,group.by = 'group',ident.1 = 'DKD',min.pct=0,logfc.threshold = 0)
CD8_marker$gene <- rownames(CD8_marker)
CD8_obj <- RunUMAP(CD8_obj,dims = 1:10)
DimPlot(CD8_obj)


CD4_obj <- subset(obj2,seurat_clusters %in% c(0,1))
CD4_obj <- RunUMAP(CD4_obj,dims = 1:15)
CD4_marker <- FindMarkers(CD4_obj,group.by = 'group',ident.1 = 'DKD',min.pct=0,logfc.threshold = 0)
CD4_marker$gene <- rownames(CD4_marker)



B_obj <- subset(obj2,seurat_clusters %in% c(3))
B_obj$group1 <- 0
B_obj$group1[B_obj$group=='DKD'] <- 1
B_marker <- FindMarkers(B_obj,group.by = 'group',ident.1 = 'DKD',min.pct=0,logfc.threshold = 0)
B_marker$gene <- rownames(B_marker)



NK_obj <- subset(obj2,seurat_clusters %in% c(6))
NK_obj$group1 <- 0
NK_obj$group1[NK_obj$group=='DKD'] <- 1
NK_marker <- FindMarkers(NK_obj,group.by = 'group',ident.1 = 'DKD',min.pct=0,logfc.threshold = 0)
NK_marker$gene <- rownames(NK_marker)


Monocyte_obj <- subset(obj2,seurat_clusters %in% c(7))
Monocyte_obj$group1 <- 0
Monocyte_obj$group1[Monocyte_obj$group=='DKD'] <- 1
Mono_marker <- FindMarkers(Monocyte_obj,group.by = 'group',ident.1 = 'DKD',min.pct=0)
Mono_marker$gene <- rownames(Mono_marker)


M2_obj <- subset(obj2,seurat_clusters %in% c(8))
M2_marker <- FindMarkers(M2_obj,group.by = 'group',ident.1 = 'DKD',min.pct=0)
M2_marker$gene <- rownames(M2_marker)

