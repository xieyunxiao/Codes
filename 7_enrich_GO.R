################################################
##GO

setwd("~/data/")
library(clusterProfiler)
library(ggplot2)
library(hrbrthemes)
library(ggsci)
library(ggpubr)
loadDb(file = "cavia.orgdb")

##########################################################################

DEGlist_UPgenes <- readRDS("DEGlist_UPgenes001.rds")
DEGlist_DOWNgenes <- readRDS("DEGlist_DOWNgenes001.rds")

upgene_number <- unlist(lapply(DEGlist_UPgenes,length))
downgene_number <- unlist(lapply(DEGlist_DOWNgenes,length))

DEGlist_UPgenes <- DEGlist_UPgenes[-which(upgene_number==0)]
DEGlist_DOWNgenes <- DEGlist_DOWNgenes[-which(downgene_number==0)]
#####################################################################33
runGO <- function(gene){ ##symbol
  if(length(gene)==0){
    return(NULL)
  }
  gene.id <- bitr(
    gene, fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = cavia.orgdb,
    drop = F
  )
  if(nrow(gene.id)==0){
    return(NULL)
  }
  gene_en <- gene.id$ENTREZID
  enrich.go <- enrichGO(gene = gene_en,  #待富集的基因列表
                        OrgDb = cavia.orgdb,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                        keyType = 'ENTREZID',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                        ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                        pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                        qvalueCutoff = 0.2,  #指定 q 值阈值（可指定 1 以输出全部）
                        readable = FALSE)
  return(enrich.go)
}
# k <- runGO(DEGlist_UPgenes$Microglia)

go_up <- lapply(DEGlist_UPgenes,runGO)
saveRDS(go_up,"enrichGOlist_UPgenes.rds")
filnames1 <- paste(names(go_up),"_GO_UP_barplot.png",sep = "")
filnames2 <- paste(names(go_up),"_GO_UP_dotplot.png",sep = "")

for (i in 1:length(go_up)) {
  if(nrow(go_up[[i]])>0){
    barplot(go_up[[i]], split="ONTOLOGY",showCategory = 6,label_format=50)+ 
      facet_grid(ONTOLOGY~.,scale="free")+
      theme(panel.grid = element_blank())+
      ggtitle(paste("Up Regulated Function in",names(go_up)[i]))+ 
      theme(plot.title = element_text(size=12,hjust=0.5))
    
    ggsave(filnames1[i],width = 8,height = 10,dpi = 600)
    
    dotplot(go_up[[i]], split="ONTOLOGY",showCategory = 6,label_format=50)+ 
      facet_grid(ONTOLOGY~.,scale="free")+
      theme(panel.grid = element_blank())+#修改主题
      theme(axis.title =element_text(size = 12, color = 'black'),
            axis.text.y =element_text(size = 12),
            legend.title=element_text(size=12))+
      scale_color_gradient(high="#FC8D62",low="#4b5cc4")+
      ggtitle(paste("Up Regulated Function in",names(go_up)[i]))+ 
      theme(plot.title = element_text(size=12,hjust=0.5))
    
    ggsave(filnames2[i],width = 8,height = 10,dpi = 600)
  }
}



go_down <- lapply(DEGlist_DOWNgenes,runGO)
saveRDS(go_down,"enrichGOlist_DOWNgenes.rds")
filnames1 <- paste(names(go_down),"_GO_DOWN_barplot.png",sep = "")
filnames2 <- paste(names(go_down),"_GO_DOWN_dotplot.png",sep = "")

for (i in 1:length(go_down)) {
  if(nrow(go_down[[i]])>0){
    barplot(go_down[[i]], split="ONTOLOGY",showCategory = 6,label_format=50)+ 
      facet_grid(ONTOLOGY~.,scale="free")+
      theme(panel.grid = element_blank())+
      ggtitle(paste("Down Regulated Function in",names(go_down)[i]))+ 
      theme(plot.title = element_text(size=12,hjust=0.5))
    
    ggsave(filnames1[i],width = 8,height = 10,dpi = 600)
    
    dotplot(go_down[[i]], split="ONTOLOGY",showCategory = 6,label_format=50)+ 
      facet_grid(ONTOLOGY~.,scale="free")+
      theme(panel.grid = element_blank())+#修改主题
      theme(axis.title =element_text(size = 12, color = 'black'),
            axis.text.y =element_text(size = 12),
            legend.title=element_text(size=12))+
      scale_color_gradient(high="#FC8D62",low="#4b5cc4")+
      ggtitle(paste("Down Regulated Function in",names(go_down)[i]))+ 
      theme(plot.title = element_text(size=12,hjust=0.5))
    
    ggsave(filnames2[i],width = 8,height = 10,dpi = 600)
  }
}

#################################################################################################3
go_down <- lapply(go_down,as.data.frame)
for (i in 1:length(go_down)) {
  if(nrow(go_down[[i]])>0){
    go_down[[i]]$CellType <- names(go_down)[i]
  }
}

go_up <- lapply(go_up,as.data.frame)
for (i in 1:length(go_up)) {
  if(nrow(go_up[[i]])>0){
    go_up[[i]]$CellType <- names(go_up)[i]
  }
}

go_up <- do.call("rbind", go_up)
go_down <- do.call("rbind", go_down)

write.table(go_up,"GO_enrich_UPgene.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(go_down,"GO_enrich_DOWNgene.txt",col.names = T,row.names = F,sep = "\t",quote = F)
############################################################################################################333
library(annotate)
df=go_up
entriz <- df$geneID
rmind <- which(entriz=="")
if(length(rmind>0)){
  df <- df[-rmind,]
}
entriz <- df$geneID
this_id <- strsplit(entriz,"/")
name <- lapply(this_id, function(x){getSYMBOL(as.character(x),data="cavia.orgdb")})
smb <- lapply(name, function(x){paste(x,collapse = "/")})
smb <- unlist(smb)

if("CellType" %in% colnames(df)){
  df <- df
}else{
  df$CellType="ALL"
}

this_df <- data.frame(CellType=df$CellType,
                      ONTOLOGY = df$ONTOLOGY,         
                      ID=df$ID,Description =df$Description,
                      p_adj=df$p.adjust,geneName=smb,geneNumber=df$Count)
write.table(this_df,"GO_enrich_UPgene.csv",col.names = T,row.names = F,quote = F,sep = ",")


#########################################################################
df=go_down
entriz <- df$geneID
rmind <- which(entriz=="")
if(length(rmind>0)){
  df <- df[-rmind,]
}
entriz <- df$geneID
this_id <- strsplit(entriz,"/")
name <- lapply(this_id, function(x){getSYMBOL(as.character(x),data="cavia.orgdb")})
smb <- lapply(name, function(x){paste(x,collapse = "/")})

smb <- unlist(smb)

if("CellType" %in% colnames(df)){
  df <- df
}else{
  df$CellType="ALL"
}
this_df <- data.frame(CellType=df$CellType,ONTOLOGY = df$ONTOLOGY, 
                      ID=df$ID,Description =df$Description,p_adj=df$p.adjust,
                      geneName=smb,geneNumber=df$Count)

write.table(this_df,"GO_enrich_DOWNgene.csv",col.names = T,row.names = F,quote = F,sep = ",")





