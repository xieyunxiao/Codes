##KEGG
setwd("~/data")
library(clusterProfiler)
library(ggplot2)
library(hrbrthemes)
library(ggsci)
library(ggpubr)
loadDB(file = "cavia.orgdb")

DEGlist_UPgenes <- readRDS("DEGlist_UPgenes001.rds")
DEGlist_DOWNgenes <- readRDS("DEGlist_DOWNgenes001.rds")

upgene_number <- unlist(lapply(DEGlist_UPgenes,length))
downgene_number <- unlist(lapply(DEGlist_DOWNgenes,length))

DEGlist_UPgenes <- DEGlist_UPgenes[-which(upgene_number==0)]
DEGlist_DOWNgenes <- DEGlist_DOWNgenes[-which(downgene_number==0)]
#####################################################################33
runKEGG <- function(gene){ ##symbol
  if(length(gene)==0){
    return(NULL)
  }
  gene.id <- bitr(
    gene, fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = cavia.orgdb
  )
  if(nrow(gene.id)==0){
    return(NULL)
  }
  
  gene_en <- gene.id$ENTREZID
  kegg <- enrichKEGG(
    gene = gene_en,  
    keyType = 'kegg',  
    organism = 'cavia.orgdb',  
    pAdjustMethod = 'fdr',  
    pvalueCutoff = 1,  
    qvalueCutoff = 1)  
  return(kegg)
}

kegg_up <- lapply(DEGlist_UPgenes,runKEGG)

saveRDS(kegg_up,"enrichKEGGlist_UPgenes_.rds")
kegg_up <- readRDS("enrichKEGGlist_UPgenes_.rds")
filnames1 <- paste(names(kegg_up),"_KEGG_UP_barplot.png",sep = "")
filnames2 <- paste(names(kegg_up),"_KEGG_UP_dotplot.png",sep = "")

for (i in 1:length(kegg_up)) {
  if(nrow(kegg_up[[i]])>0){
    barplot(kegg_up[[i]], showCategory =10)+ggtitle(paste("Up Regulated Pathway in",names(kegg_up)[i]))+ 
      theme(plot.title = element_text(size=12,hjust=0.5))
    ggsave(filnames1[i],width = 8,height = 6,dpi = 600)
    
    dotplot(kegg_up[[i]],x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
            showCategory = 10) +ggtitle(paste("Up Regulated Pathway in ",names(kegg_up)[i]))+
      theme(plot.title = element_text(size=12,hjust=0.5))
    ggsave(filnames2[i],width = 8,height = 6,dpi = 600)
  }
}

kegg_down <- lapply(DEGlist_DOWNgenes,runKEGG)

saveRDS(kegg_down,"enrichKEGGlist_DOWNgenes.rds")
kegg_down <- readRDS("enrichKEGGlist_DOWNgenes.rds")
filnames1 <- paste(names(kegg_down),"_KEGG_DOWN_barplot.png",sep = "")
filnames2 <- paste(names(kegg_down),"_KEGG_DOWN_dotplot.png",sep = "")

for (i in 1:length(kegg_down)) {
  if(is.null(kegg_down[[i]]) || nrow(kegg_down[[i]])==0){
    next
  }
  barplot(kegg_down[[i]], showCategory =10)+ggtitle(paste("Down Regulated Pathway in",names(kegg_down)[i]))+ 
    theme(plot.title = element_text(size=12,hjust=0.5))
  ggsave(filnames1[i],width = 8,height = 6,dpi = 600)
  
  dotplot(kegg_down[[i]],x = "GeneRatio", color = "p.adjust", size = "Count", #默认参数
          showCategory = 10) +ggtitle(paste("Down Regulated Pathway in ",names(kegg_down)[i]))+
    theme(plot.title = element_text(size=12,hjust=0.5))
  ggsave(filnames2[i],width = 8,height = 6,dpi = 600)
}

############################################################################################################33
kegg_down <- lapply(kegg_down,as.data.frame)
for (i in 1:length(kegg_down)) {
  if(nrow(kegg_down[[i]])>0){
    kegg_down[[i]]$CellType <- names(kegg_down)[i]
  }
}

kegg_up <- lapply(kegg_up,as.data.frame)
for (i in 1:length(kegg_up)) {
  if(nrow(kegg_up[[i]])>0){
    kegg_up[[i]]$CellType <- names(kegg_up)[i]
  }
}

kegg_up <- do.call("rbind", kegg_up)
kegg_down <- do.call("rbind", kegg_down)

write.table(kegg_up,"KEGG_enrich_UPgene.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(kegg_down,"KEGG_enrich_DOWNgene.txt",col.names = T,row.names = F,sep = "\t",quote = F)
#######################################################################################################33

############################################################################################################333
library(annotate)
df=kegg_up
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
                      ID=df$ID,Description =df$Description,
                      p_adj=df$p.adjust,geneName=smb,geneNumber=df$Count)

write.table(this_df,"KEGG_enrich_UPgene.csv",col.names = T,row.names = F,quote = F,sep = ",")

#########################################################################
df=kegg_down
entriz <- df$geneID
rmind <- which(entriz=="")
if(length(rmind>0)){
  df <- df[-rmind,]
}
entriz <- df$geneID
this_id <- strsplit(entriz,"/")
name <- lapply(this_id, function(x){getSYMBOL(as.character(x),data="org.Hs.eg.db")})
smb <- lapply(name, function(x){paste(x,collapse = "/")})

smb <- unlist(smb)

if("CellType" %in% colnames(df)){
  df <- df
}else{
  df$CellType="ALL"
}
this_df <- data.frame(CellType=df$CellType,ID=df$ID,
                      Description =df$Description,
                      p_adj=df$p.adjust,geneName=smb,
                      geneNumber=df$Count)

write.table(this_df,"KEGG_enrich_DOWNgene.csv",col.names = T,row.names = F,quote = F,sep = ",")

###########################################################################3
