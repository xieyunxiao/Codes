##plot
##############################################################################################3
## cell ratio plot
sam2ct <- pbmc@meta.data[,c("sampleinf","Cell_Type")]
sam2ct_stat <- as.matrix(table(sam2ct))
sam2ct_ratio <- sam2ct_stat/rowSums(sam2ct_stat)
ratio <- as.data.frame(sam2ct_ratio)
colnames(ratio) <- c("Group","Cell_Type","Cell_Ratio")
lev <- levels(ratio$Cell_Type)
ratio$Cell_Type <- as.character(ratio$Cell_Type)

ratio$Cell_Type <- factor(ratio$Cell_Type,levels = lev,labels = lev)
LIM <- ratio[which(ratio$Group=="LIM"),]
NC <- ratio[which(ratio$Group=="NC"),]
a <- rev(LIM$Cell_Ratio)
b <- rev(LIM$Cell_Type)
lab <- paste0(round(a * 100, 2), "%")
ang <- rep(0,length(lab))

ggplot(LIM,aes(x = 1, y = Cell_Ratio, fill = Cell_Type)) +
  geom_col(colour = "white")+coord_cartesian (xlim = c (0, 5))+
  geom_text(aes(y = as.numeric(a/2 + c(0, cumsum(a)[-length(a)])),
                x = 1.8, label = lab),
            hjust=0.5,vjust=0.5,size=5,angle=ang)+
  coord_polar(theta = "y", start = 1.65) +ggtitle("LIM") +
  xlim(c(-1, 2)) +
  theme(plot.title = element_text(hjust = 0.5,size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  ) 
ggsave("LIM_CellRate_label.png",width = 7,height = 7,dpi = 600)

########################################################################
a <- rev(NC$Cell_Ratio)
b <- rev(NC$Cell_Type)
lab <- paste0(round(a * 100, 2), "%")
ang <- rep(0,length(lab))
ggplot(NC,aes(x = 1, y = Cell_Ratio, fill = Cell_Type)) +
  geom_col(colour = "white")+coord_cartesian (xlim = c (0, 5))+
  geom_text(aes(y = as.numeric(a/2 + c(0, cumsum(a)[-length(a)])),
                x = 1.8, label = lab),
            hjust=0.5,vjust=0.5,size=5,angle=ang)+
  coord_polar(theta = "y", start = 1.65) +ggtitle("NC") +
  xlim(c(-1, 2)) +
  theme(plot.title = element_text(hjust = 0.5,size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
  ) 
ggsave("NC_CellRate_label.png",width = 7,height = 7,dpi = 600)


###################################################################################################################33
DimPlot(pbmc,cols = rep("grey",11))+NoLegend()+ggtitle("All Cells")+theme(plot.title = element_text(hjust = 0.5))
ggsave("All_Cells.png",width = 5,height = 5,dpi = 600)
DimPlot(pbmc,group.by = "RNA_snn_res.0.2",cols = rainbow(19))+NoLegend()+ggtitle("Pre-Cluster")+theme(plot.title = element_text(hjust = 0.5))
ggsave("Pre_cluster.png",width = 5,height = 5,dpi = 600)
DimPlot(pbmc,label = F)+NoLegend()+ggtitle("Cell Annotation")+theme(plot.title = element_text(hjust = 0.5))
ggsave("Annotation.png",width = 5,height = 5,dpi = 600)

DimPlot(pbmc,label = T,label.size = 5,group.by="main_class")
  ggtitle("Cell Main Class")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("CellMainClass2.png",width = 7,height = 5,dpi = 600)

colour <- c("#F8766D", "#E7851E", "#D09400", "#B2A100", "#89AC00", "#45B500", "#00BC51", "#00C087", 
            "#00C0B2", "#00BCD6", "#00B3F2", "#29A3FF", "#9C8DFF", "#D277FF", "#F166E8", "#EA9999", "#A9689A")
                        
DimPlot(pbmc,label = T,repel = T,label.size=5,cols = colour)+
  ggtitle("Cell Type")+theme(plot.title = element_text(hjust = 0.5))
ggsave("CellType.png",width = 6,height = 5,dpi = 600)

DimPlot(pbmc,label = T,label.size=5,split.by="stim",group.by="main_class")+
  NoLegend()+ggtitle("Cell Main Class")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave("CellMainClass_sample.png",width = 9,height = 5,dpi = 600)

DimPlot(pbmc,label = T,repel = T,label.size=5,split.by="stim",cols = colour)+ggtitle("Cell Type")+theme(plot.title = element_text(hjust = 0.5))
ggsave("CellType_sample.png",width = 10,height = 5,dpi = 600)
