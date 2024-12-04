##Cell annotation
setwd("~/data")
library(Seurat)

vcall <- readRDS("VC_combined_snsc.rds")
Idents(vcall) <- vcall$RNA_snn_res.0.3
vcall.markers <- FindAllMarkers(vcall, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(vcall.markers,"vcall_markers_cluster.rds")

Idents(vcall) <- vcall$RNA_snn_res.0.3

FeaturePlot(vcall, features = toupper(c("Cux2","Ccbe1","Mdga1","Stard8"))) ##L2L3
FeaturePlot(vcall, features = toupper(c("Whrn","Rorb"))) ##L4
FeaturePlot(vcall, features = toupper(c("Rorb","Deptor","Foxo1","Ptprm")))##L5IP
FeaturePlot(vcall, features = toupper(c("Nxph1","Tshz2","Trhr","Slc17a8"))) # L5NP
FeaturePlot(vcall, features = toupper(c("Bcl6","Erg","Reln","Zfp804b"))) #L5PT
FeaturePlot(vcall, features = toupper(c("Foxp2","Syt6","Cdh9")))#L6CT
FeaturePlot(vcall, features = toupper(c("Zfp804b","Cdh9")))#L6IT
FeaturePlot(vcall, features = toupper(c("Ctgf","Inpp4b","Svil","Foxp2")))#L6b
FeaturePlot(vcall, features = c("NDNF","VIP","SST","PVALB"))  ##GABA
FeaturePlot(vcall, features = c("HSD11B1","TCERG1L","PDE1C","BATF3"))  ##GABA PVALB
FeaturePlot(vcall, features = toupper(c("GFAP","Aldh1l1","Atp1b2","Camsap1"))) ##Astrocyte
FeaturePlot(vcall, features = c("OLIG2","CNP","AF2H","UGTB","MBP","PTGDS")) ##Oligodendrocyte
FeaturePlot(vcall, features = c("PCDH15","LHFPL3"))##OPC
FeaturePlot(vcall, features = c("IBA1","AIF1","CD68","CX3CR1","DOCK2")) # Microglia
FeaturePlot(vcall, features = c("COL1A2"),max.cutoff = 2) #Vascular Cell
FeaturePlot(vcall, features = c("CD31","PECAM1","CD34","CD309","CD133"))#Endotheial
FeaturePlot(vcall, features = c("DOCK2","CD68"))#Macropage
FeaturePlot(vcall, features = c("MRC1","CD163","STAB1","CTSC","CPM"))#Macropage M2
FeaturePlot(vcall, features = c("SLC17A7","RORB","GAD1","LHX6","NR2F2"))  ##Neurons
FeaturePlot(vcall, features = c("RGS5","ABCC9"))  ##smc
FeaturePlot(vcall, features = c("CXCL8","S100A9","G0S2","LCP1")) ##Monocytes & Neutrophils
FeaturePlot(vcall, features = c("NKG7","CD3E","SATB1","SKAP1")) ##T-cells & NK-cells
FeaturePlot(vcall, features = c("TOP2A","KIF11","UBE2C","CENPF")) ##NPC
FeaturePlot(vcall, features = c("CD44")) ##stem cell
FeaturePlot(vcall, features = c("CD44","S100A4")) ##mesenchymal stem cell
FeaturePlot(vcall, features = c("GAD1", "GAD2", "DLX2")) ##

################################################################################################################################333
cl <- 0:20
ct <- c("L2 & L3","Microglia","L6","Oligodendrocyte","Astrocyte","L6","L6",
        "SMC","L5","L5","OPC","L5","Microglia","Macrophage","Oligodendrocyte","L5",
        "Astrocyte","SMC","OPC","Monocytes & Neutrophils","Astrocyte")
vcall$Cell_Type <- factor(vcall$RNA_snn_res.0.3,levels = cl,labels = ct)
vcall$Cell_Type <- as.character(vcall$Cell_Type)
vcall$Cell_Type[which(vcall$RNA_snn_res.0.7 %in% c(4,13))] <- "L4"

ct <- c("L2 & L3","L4","L5","L6","Astrocyte","Microglia","OPC","Oligodendrocyte",
        "SMC","Macrophage","Monocytes & Neutrophils")
vcall$Cell_Type <- factor(vcall$Cell_Type,levels = ct,labels = ct)

Idents(vcall) <- vcall$Cell_Type

DimPlot(vcall, reduction = "umap", label = F, pt.size = 0.01,split.by = "sampleinf",group.by = "Cell_Type")+ggtitle("Sample Info")

DimPlot(vcall, reduction = "umap",  pt.size = 0.01,group.by = "Cell_Type")+ggtitle("Cell Type")

DimPlot(vcall, reduction = "umap", label = T, pt.size = 0.01,group.by = "Cell_Type",repel = T)+ggtitle("Cell Type")
saveRDS(vcall,"VC_combined_snsc.rds")

##################################################################################################33
##############################################################################################3
#tracks plot

library(vctrs)
library(scRNAtoolVis)
# load test data

pbmc$Cell_Type <- as.character(pbmc$Cell_Type)

ct <- levels(pbmc$Cell_Type)

Idents(pbmc) <- pbmc$Cell_Type
table(pbmc$Cell_Type)
table(pbmc$main_class)

markerslist <- vcall.markers
markerslist <- split(markerslist,markerslist$cluster)
names(markerslist)
cellMarkers <- list(SST=c("SST","GRIN3A","RELN"),
                    VIP=c("VIP","PROX1","CRH","GRPR"),
                    PVALB=c("PVALB","TMEM132C"),
                    L2_L3=c("CUX2","CCBE1","STARD8"),
                    L4=c("WHRN","RORB","CTXN3"),
                    L5=c("TSHZ2"),
                    L5_Pde1c = c("PDE1C"),
                    L6=c("FOXP2","CRYM","SYT6"),
                    Astro=toupper(c("GFAP","Aldh1l1","Atp1b2")),
                    Microglia=c("AIF1","CX3CR1","DOCK2"),
                    OPC=c("PCDH15","LHFPL3"),
                    Oligo=c("CNP","MBP","PTGDS"),
                    SMC=c("RGS5"),
                    Macrophage=c("MRC1","CD163","STAB1"),
                    Mon_Neutr=c("CXCL8","S100A9","G0S2","LCP1")
)


markerslist2 <- markerslist
for (i in 1:length(cellMarkers)) {
  m<- markerslist[[i]][cellMarkers[[i]],]
  m <- m[!is.na(m[,1]),]
  markerslist2[[i]] <- m
}
dup <- unlist(lapply(markerslist2,nrow))
markers2 <- do.call("rbind", markerslist2)
markers2 <- markers2[which(markers2$avg_log2FC>0),]

saveRDS(markers2,"markersInCellMarkers2.rds")

markers2 <- markers2[!duplicated(markers2$gene),]

# tracksplot color
ob2 <- subset(x = pbmc, downsample = 150)
ob2$Cell_Type <- as.character(ob2$Cell_Type)
table(ob2$Cell_Type)

ob2$Cell_Type <- factor(ob2$Cell_Type,levels = names(markerslist),labels = names(markerslist))
Idents(ob2) <- ob2$Cell_Type
table(ob2$Cell_Type)
features <- unlist(cellMarkers)
features <- as.vector(features)
markers2$gene <- as.character(markers2$gene)
markers2$cluster <- as.character(markers2$cluster)
table(markers2$cluster)
markers2$cluster <- factor(markers2$cluster,levels = names(markerslist),labels = names(markerslist))


p <- tracksPlot(object = ob2,
           plot.type = "track",
           genes = markers2,
           strip_nested_params = list(background_x = ggh4x::elem_list_rect(fill = rainbow(15))),
           theme_params = list(panel.spacing.x = unit(0,"mm")),
           # facet_nested_params = list(space = "free"),
           cell.order = names(markerslist),
           gene.order = as.vector(features))+NoLegend()+
  scale_fill_manual(values = rep("grey",length(features)))


pdf("CellMarkers.pdf",width = 23,height = 10)
p+NoLegend()+
  scale_fill_manual(values = rep(rep(rainbow(15),dup),length(features)))
dev.off()
