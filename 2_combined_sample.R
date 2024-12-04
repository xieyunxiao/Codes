setwd("~/data")
library(Seurat)
library(hrbrthemes)
library(ggsci)
library(ggpubr)
setwd("~/VC/")
vcsn <- readRDS("Single_Nuclei_RNA_Seq//sn_combined.rds")
vcsc <- readRDS("Single_Cell_RNA_Seq//sc_combined.rds")

table(vcsc$orig.ident)
table(vcsn$orig.ident)

vcall<- merge(vcsn, y = vcsc)
Idents(vcall) <- vcall$CellType

rm(vcsn)
rm(vcsc)
#NormalizeData?
vcall <- NormalizeData(vcall)

##FindVariableFeatures
vcall <- FindVariableFeatures(vcall, selection.method = "vst", nfeatures = 2000)

#ScaleData, vcsn[["RNA"]]@scale.data
all.genes <- rownames(vcall)
write.table(all.genes,"allgene_snsc.txt",row.names = F,col.names = F,quote = F)
vcall <- ScaleData(vcall, features = all.genes)

#PCA
vcall <- RunPCA(vcall,features = VariableFeatures(object = vcall),npcs = 50)


##PC number
vcall <- JackStraw(vcall, dims = 50,num.replicate = 100)
vcall <- ScoreJackStraw(vcall, dims = 1:50)
JackStrawPlot(vcall, dims = 1:15)
y <- ElbowPlot(vcall,ndims = 30)
y
# y$data

##FindClusters
vcall <- FindNeighbors(vcall, dims = 1:30)
for (i in seq(0.1,1,by=0.1)) {
  vcall <- FindClusters(vcall, resolution = i)
}
# head(vcall@meta.data)
# vcall$seurat_clusters <- vcall$RNA_snn_res.0.3
Idents(vcall) <- vcall$CellType
head(Idents(vcall),5)

##UMAP/tSNE
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
vcall <- RunUMAP(vcall, dims = 1:30)
vcall <- RunTSNE(vcall, dims = 1:30)

vcall$batch <- "sn"
vcall$batch[which(vcall$orig.ident %in% c("VCscNC4W","VCscLIM4W"))] <- "sc"

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(vcall, reduction = "umap")
DimPlot(vcall, reduction = "umap",group.by = "orig.ident")
DimPlot(vcall, reduction = "umap",split.by = "sampleinf")
DimPlot(vcall, reduction = "tsne")
DimPlot(vcall, reduction = "tsne",group.by="sampleinf")
DimPlot(vcall, reduction = "tsne",split.by="sampleinf")

DimPlot(vcall, reduction = "umap",group.by = "sampleinf",label = T) + ggtitle("Before remove batch effect")
DimPlot(vcall, reduction = "tsne",group.by = "sampleinf",label = T) + ggtitle("Before remove batch effect")

saveRDS(vcall,"VC_combined_snsc.rds")
