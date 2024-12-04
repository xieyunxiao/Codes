##DEGs suitable for both VC and Retina
ctps <- levels(pbmc$Cell_Type)
markers_all <- FindMarkers(pbmc, ident.1 = "LIM",ident.2="NC",  min.pct = 0.25,group.by = "sampleinf")
markers_all$cluster <- "ALL"
markers_all$gene <- rownames(markers_all)

markerslist <-list()
for (i in 1:length(ctps)) {
  subvcsn <- subset(pbmc,subset = Cell_Type == ctps[i])
  markers <- FindMarkers(subvcsn, ident.1 = "LIM",ident.2="NC",  min.pct = 0.25,group.by = "sampleinf")
  markers$cluster <- ctps[i]
  markers$gene <- rownames(markers)
  markerslist[[ctps[i]]] <- markers
}

markerslist[["All"]] <- markers_all
markerslist2 <- do.call("rbind", markerslist)
tp <- unique(markerslist2$cluster)
names(markerslist) <- tp
saveRDS(markerslist,"LIM2NCdiffGeneList.rds")

markerslist2 <- do.call("rbind", markerslist)
saveRDS(markerslist2,"markerslist2.rds")