library(Seurat)
library(hrbrthemes)
library(ggsci)
library(ggpubr)

### VC Single Nuclei RNA-Seq
setwd("~/VC/Single_Nuclei_RNA_Seq/")
NC <- Read10X(data.dir = "NC_4w_sn/")
LIM <- Read10X(data.dir = "LIM_4w_sn/")
# Initialize the Seurat object with the raw (non-normalized data).
NC <- CreateSeuratObject(counts = NC, project = "NC4Wsn", min.cells = 3, min.features = 200)
LIM <- CreateSeuratObject(counts = LIM, project = "LIM4Wsn", min.cells = 3, min.features = 200)
NC$stim  <- "NC"
NC$sampleinf  <- "NC"
LIM$stim  <- "LIM"
LIM$sampleinf  <- "LIM"
pbmc<- merge(NC, y = LIM)
Idents(pbmc) <- pbmc$orig.ident

##PercentageFeatureSet() QC 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))

plot2 <- VlnPlot(pbmc, features = c( "nCount_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))
  
plot3 <- VlnPlot(pbmc, features = c("percent.mt"),raster = F)+NoLegend()
plot4 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F)+
  scale_x_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))+
  theme(axis.text.x = element_text(angle=30))

# plot1+plot2+plot3
ggarrange(plot1, plot2, plot3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
plot4

#subset
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & percent.mt < 5 & nCount_RNA > 2000)

plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))

plot2 <- VlnPlot(pbmc, features = c( "nCount_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))

plot3 <- VlnPlot(pbmc, features = c("percent.mt"),raster = F)+NoLegend()
plot4 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F)+
  scale_x_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))+
  theme(axis.text.x = element_text(angle=30))

# plot1+plot2+plot3
ggarrange(plot1, plot2, plot3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
plot4


#NormalizeData?
pbmc <- NormalizeData(pbmc)


##FindVariableFeatures
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#ScaleData. saved in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc),npcs = 50)

DimPlot(pbmc)
DimHeatmap(pbmc)
# 
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


##PCA PC number
pbmc <- JackStraw(pbmc, dims = 50,num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
JackStrawPlot(pbmc, dims = 1:15)
y <- ElbowPlot(pbmc,ndims = 30)
y
# y$data

##FindClusters
pbmc <- FindNeighbors(pbmc, dims = 1:20)
for (i in seq(0.1,1,by=0.1)) {
  pbmc <- FindClusters(pbmc, resolution = i)
}
head(pbmc@meta.data)
pbmc$seurat_clusters <- pbmc$RNA_snn_res.1
Idents(pbmc) <- pbmc$seurat_clusters
head(Idents(pbmc),5)

##UMAP/tSNE
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages ='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",group.by = "sampleinf")
DimPlot(pbmc, reduction = "umap",split.by = "sampleinf")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",group.by="sampleinf")
DimPlot(pbmc, reduction = "tsne",split.by="sampleinf")
saveRDS(pbmc, file = "sn_combined.rds")

###########################################################################################################################################################
##########################################################################################################################################################
### VC Single Cell RNA-Seq
setwd("~/VC/Single_Cell_RNA_Seq/")
NC <- Read10X(data.dir = "NC_4w_sc/")
LIM <- Read10X(data.dir = "LIM_4w_sc/")
# Initialize the Seurat object with the raw (non-normalized data).
NC <- CreateSeuratObject(counts = NC, project = "NC4Wsc", min.cells = 3, min.features = 200)
LIM <- CreateSeuratObject(counts = LIM, project = "LIM4Wsc", min.cells = 3, min.features = 200)
NC$stim  <- "NC"
NC$sampleinf  <- "NC"
LIM$stim  <- "LIM"
LIM$sampleinf  <- "LIM"
pbmc<- merge(NC, y = LIM)
Idents(pbmc) <- pbmc$orig.ident

##PercentageFeatureSet() QC 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))

plot2 <- VlnPlot(pbmc, features = c( "nCount_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))
  
plot3 <- VlnPlot(pbmc, features = c("percent.mt"),raster = F)+NoLegend()
plot4 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F)+
  scale_x_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))+
  theme(axis.text.x = element_text(angle=30))

# plot1+plot2+plot3
ggarrange(plot1, plot2, plot3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
plot4

#subset
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & percent.mt < 5 & nCount_RNA > 2000)

plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))

plot2 <- VlnPlot(pbmc, features = c( "nCount_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))

plot3 <- VlnPlot(pbmc, features = c("percent.mt"),raster = F)+NoLegend()
plot4 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F)+
  scale_x_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))+
  theme(axis.text.x = element_text(angle=30))

# plot1+plot2+plot3
ggarrange(plot1, plot2, plot3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
plot4


#NormalizeData?
pbmc <- NormalizeData(pbmc)


##FindVariableFeatures
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#ScaleData. saved in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc),npcs = 50)

DimPlot(pbmc)
DimHeatmap(pbmc)
# 
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


##PCA PC number
pbmc <- JackStraw(pbmc, dims = 50,num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
JackStrawPlot(pbmc, dims = 1:15)
y <- ElbowPlot(pbmc,ndims = 30)
y
# y$data

##FindClusters
pbmc <- FindNeighbors(pbmc, dims = 1:20)
for (i in seq(0.1,1,by=0.1)) {
  pbmc <- FindClusters(pbmc, resolution = i)
}
head(pbmc@meta.data)
pbmc$seurat_clusters <- pbmc$RNA_snn_res.1
Idents(pbmc) <- pbmc$seurat_clusters
head(Idents(pbmc),5)

##UMAP/tSNE
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages ='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",group.by = "sampleinf")
DimPlot(pbmc, reduction = "umap",split.by = "sampleinf")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",group.by="sampleinf")
DimPlot(pbmc, reduction = "tsne",split.by="sampleinf")
saveRDS(pbmc, file = "sc_combined.rds")


###########################################################################################################################################################
##########################################################################################################################################################
### Retina Single Cell RNA-Seq
setwd("~/Retina/Single_Cell_RNA_Seq/")
NC <- Read10X(data.dir = "NC_4w_R/")
LIM <- Read10X(data.dir = "LIM_4w_R/")
# Initialize the Seurat object with the raw (non-normalized data).
NC <- CreateSeuratObject(counts = NC, project = "NC4W", min.cells = 3, min.features = 200)
LIM <- CreateSeuratObject(counts = LIM, project = "LIM4W", min.cells = 3, min.features = 200)
NC$stim  <- "NC"
NC$sampleinf  <- "NC"
LIM$stim  <- "LIM"
LIM$sampleinf  <- "LIM"
pbmc<- merge(NC, y = LIM)
Idents(pbmc) <- pbmc$orig.ident

##PercentageFeatureSet() QC 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))

plot2 <- VlnPlot(pbmc, features = c( "nCount_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))
  
plot3 <- VlnPlot(pbmc, features = c("percent.mt"),raster = F)+NoLegend()
plot4 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F)+
  scale_x_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))+
  theme(axis.text.x = element_text(angle=30))

# plot1+plot2+plot3
ggarrange(plot1, plot2, plot3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
plot4

#subset
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & percent.mt < 5 & nCount_RNA > 2000)

plot1 <- VlnPlot(pbmc, features = c("nFeature_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))

plot2 <- VlnPlot(pbmc, features = c( "nCount_RNA"),raster = F)+NoLegend()+
  scale_y_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))

plot3 <- VlnPlot(pbmc, features = c("percent.mt"),raster = F)+NoLegend()
plot4 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F)+
  scale_x_continuous(breaks = seq(0,max(pbmc$nCount_RNA),by=5000))+
  scale_y_continuous(breaks = seq(0,max(pbmc$nFeature_RNA),by=500))+
  theme(axis.text.x = element_text(angle=30))

# plot1+plot2+plot3
ggarrange(plot1, plot2, plot3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
plot4


#NormalizeData?
pbmc <- NormalizeData(pbmc)


##FindVariableFeatures
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

#ScaleData. saved in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc),npcs = 50)

DimPlot(pbmc)
DimHeatmap(pbmc)
# 
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


##PCA PC number
pbmc <- JackStraw(pbmc, dims = 50,num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
JackStrawPlot(pbmc, dims = 1:15)
y <- ElbowPlot(pbmc,ndims = 30)
y
# y$data

##FindClusters
pbmc <- FindNeighbors(pbmc, dims = 1:20)
for (i in seq(0.1,1,by=0.1)) {
  pbmc <- FindClusters(pbmc, resolution = i)
}
head(pbmc@meta.data)
pbmc$seurat_clusters <- pbmc$RNA_snn_res.1
Idents(pbmc) <- pbmc$seurat_clusters
head(Idents(pbmc),5)

##UMAP/tSNE
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages ='umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",group.by = "sampleinf")
DimPlot(pbmc, reduction = "umap",split.by = "sampleinf")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",group.by="sampleinf")
DimPlot(pbmc, reduction = "tsne",split.by="sampleinf")
DimPlot(pbmc, reduction = "umap",group.by = "sampleinf",label = T) + ggtitle("Before remove batch effect")
DimPlot(pbmc, reduction = "tsne",group.by = "sampleinf",label = T) + ggtitle("Before remove batch effect")
saveRDS(pbmc, file = "Retina_combined.rds")