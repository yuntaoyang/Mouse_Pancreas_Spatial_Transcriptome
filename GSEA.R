#---- Load Packages ------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

#---- Set Directory ------------------------------------------------------------
out_dir <- './out/'

#---- Pathway Analysis ---------------------------------------------------------
# preprocessing
aggr_AAP_vs_YAP <- read.csv(paste0(out_dir,'aggr_AAP_vs_YAP.csv'), row.names = 'X')
aggr_AAP_vs_YAP$gene <- row.names(aggr_AAP_vs_YAP)
ids <- bitr(aggr_AAP_vs_YAP$gene,'SYMBOL','ENTREZID','org.Mm.eg.db')
data <- merge(aggr_AAP_vs_YAP,ids,by.x='gene',by.y='SYMBOL')
data <- data[order(data$avg_log2FC,decreasing = T),]
# GO-BP
list_symbol <- as.numeric(data$avg_log2FC)
names(list_symbol) <- data$gene
gseaGO <- gseGO(list_symbol,OrgDb=org.Mm.eg.db,
                ont='BP',keyType="SYMBOL",pAdjustMethod="BH",
                minGSSize=10,maxGSSize=500,eps=0,
                pvalueCutoff=1,verbose=FALSE,by="fgsea")
write.csv(gseaGO,paste0(out_dir,'aggr_AAP_vs_YAP_GOBP.csv'),row.names = FALSE)
# KEGG
list_id <- as.numeric(data$avg_log2FC)
names(list_id) <- data$ENTREZID
egseKEGG <- gseKEGG(list_id,organism='mmu',keyType="kegg",
                    minGSSize=10,maxGSSize=500,eps=0,
                    pvalueCutoff=1,pAdjustMethod="BH",
                    verbose=FALSE,by="fgsea")
egseKEGG <- as.data.frame(egseKEGG)
egseKEGG$core_enrichment <- as.character(egseKEGG$core_enrichment)
gene_name = c()
for(i in egseKEGG$core_enrichment)
{
  gene_list = c()
  for (j in as.list(strsplit(i,'/'))[[1]])
  {
    gene_list <- c(gene_list,paste(data[data$ENTREZID==j,]$gene,collapse='/'))
  }
  gene_name <- c(gene_name,paste(gene_list,collapse='/'))
}
egseKEGG <- subset(egseKEGG, select = -c(core_enrichment))
egseKEGG['core_enrichment'] <- gene_name
write.csv(egseKEGG,paste0(out_dir,'aggr_AAP_vs_YAP_KEGG.csv'),row.names = FALSE)
