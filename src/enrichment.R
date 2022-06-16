library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(tibble)
library(ReactomePA)

deg <- read.delim("data.txt", header = T,
                  stringsAsFactors = F)
deg <- deg[deg$qval < 0.05,]

deg$entrez <- mapIds(x = org.Hs.eg.db,
                      keys = deg$geneName,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")

p1_geneList <- as.data.frame(cbind(deg$entrez,deg$log2FoldChange)) %>% deframe()

#--------------------------------------------
#GO
egoBP <- enrichGO(gene = as.list(deg$entrez), 
                  OrgDb    = org.Hs.eg.db,
                  ont      = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable = TRUE)

edox1 <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')

p1_map <- cnetplot(edox1, foldChange=p1_geneList, node_label = "all",
                    colorEdge = TRUE, categorySize="geneNum",
                    cex_category = 0.8, 
                    cex_label_gene = 0.8, 
                    cex_label_category = 0.8,
                    showCategory = 5)

tiff("GO_enrichment.tiff", width = 5000, height = 4000, res = 390) 
p1_map
dev.off()

#--------------------------------------------
#KEGG
kk <- enrichKEGG(gene = deg$entrez,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

KK2 <- setReadable(kk, OrgDb  = org.Hs.eg.db, keyType="ENTREZID")

p2_map <- cnetplot(KK2,foldChange=p1_geneList, node_label = "all",
                    colorEdge = TRUE, categorySize="geneNum",
                    cex_category = 0.8, 
                    cex_label_gene = 0.9, 
                    cex_label_category = 0.7,
                    showCategory = 5)

tiff("KEEG_Erichment.tiff", width = 5000, height = 4000, res = 390) 
p2_map
dev.off()

#-------------------------------------------------------
#Reactome
rec <- enrichPathway(gene=deg$entrez, pvalueCutoff = 0.05, 
                      readable=TRUE)

p3_map <- cnetplot(rec,foldChange=p1_geneList, node_label = "all",
                    colorEdge = TRUE, categorySize="geneNum",
                    cex_category = 0.8, 
                    cex_label_gene = 0.9, 
                    cex_label_category = 0.7,
                    showCategory = 5)
                    
tiff("Reactome_enrichment.tiff", width = 5000, height = 4000, res = 390) 
p3_map
dev.off()
