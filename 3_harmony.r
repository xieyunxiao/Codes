##harmony
setwd("~/data")
library(Seurat)
library(hrbrthemes)
library(ggsci)
library(ggpubr)

library(harmony)
vcall <- readRDS("VC_combined_snsc.rds")
vcall <- RunHarmony(vcall,group.by.vars = "orig.ident",reduction.save = "harmony",max.iter.harmony = 30)
DimPlot(vcall, reduction = "harmony")

vcall <- FindNeighbors(vcall,reduction = "harmony", dims = 1:20)
for (i in seq(0.1,1,by=0.1)) {
  vcall <- FindClusters(vcall, resolution = i)
}

vcall <- RunUMAP(vcall, dims = 1:20,reduction = "harmony")
vcall <- RunTSNE(vcall, dims = 1:20,reduction = "harmony")

vcall$UMAP_1 <- vcall@reductions$umap@cell.embeddings[,1]
vcall$UMAP_2 <- vcall@reductions$umap@cell.embeddings[,2]
vcall$tSNE_1 <- vcall@reductions$tsne@cell.embeddings[,1]
vcall$tSNE_2 <- vcall@reductions$tsne@cell.embeddings[,2]

head(vcall@meta.data)
DimPlot(vcall, reduction = "umap",group.by = "batch")
DimPlot(vcall, reduction = "umap",split.by = "batch")


DimPlot(vcall, reduction = "umap",group.by = "orig.ident",label = T) + ggtitle("After remove batch effect")
DimPlot(vcall, reduction = "tsne",group.by = "orig.ident",label = T) + ggtitle("After remove batch effect")
saveRDS(vcall, file = "VC_combined_snsc.rds")

##############################################################################################################3
pbmc <- readRDS("Retina_combined.rds")
pbmc <- RunHarmony(pbmc,group.by.vars = "orig.ident",reduction.save = "harmony",max.iter.harmony = 30)
DimPlot(pbmc, reduction = "harmony")

pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = 1:20)
for (i in seq(0.1,1,by=0.1)) {
  pbmc <- FindClusters(pbmc, resolution = i)
}

pbmc <- RunUMAP(pbmc, dims = 1:20,reduction = "harmony")
pbmc <- RunTSNE(pbmc, dims = 1:20,reduction = "harmony")

pbmc$UMAP_1 <- pbmc@reductions$umap@cell.embeddings[,1]
pbmc$UMAP_2 <- pbmc@reductions$umap@cell.embeddings[,2]
pbmc$tSNE_1 <- pbmc@reductions$tsne@cell.embeddings[,1]
pbmc$tSNE_2 <- pbmc@reductions$tsne@cell.embeddings[,2]


DimPlot(pbmc, reduction = "umap",group.by = "orig.ident",label = T) + ggtitle("After remove batch effect")
DimPlot(pbmc, reduction = "tsne",group.by = "orig.ident",label = T) + ggtitle("After remove batch effect")
saveRDS(pbmc, file = "Retina_combined.rds")
