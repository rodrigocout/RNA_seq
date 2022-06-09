if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("DEGreport")
#install.packages("/Users/rodrigoalmeida/DEGreport_1.26.0.tgz", repos = NULL, type="source")
#devtools::install_git("https://git@git.bioconductor.org/packages/DEGreport")
#install.packages("devtools")
#devtools::install_github("lpantano/DEGreport")
#install.packages("Cairo")

library("DEGreport")
data(humanGender)

library(DESeq2)
idx <- c(1:10, 75:85)
dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
                              colData(humanGender)[idx,], design=~group)
dds <- DESeq(dds)
res <- results(dds)

counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))
degCheckFactors(counts[, 1:6])

degQC(counts, design[["group"]], pvalue = res[["pvalue"]])

resCov <- degCovariates(log2(counts(dds)+0.5),
                        colData(dds))

cor <- degCorCov(colData(dds))

data(geneInfo)
degSignature(humanGender, geneInfo, group = "group")

ma = assay(rlog(dds))[row.names(res)[1:100],]
res <- degPatterns(ma, design, time = "group")

filter_count <- degFilter(counts(dds),
                          design, "group",
                          min=1, minreads = 50)
cat("gene in final count matrix", nrow(filter_count))

library(ComplexHeatmap)
th <- HeatmapAnnotation(df = colData(dds),
                        col = degColors(colData(dds), TRUE))
Heatmap(log2(counts(dds) + 0.5)[1:10,],
        top_annotation = th)

library(pheatmap)
pheatmap(log2(counts(dds) + 0.5)[1:10,], 
         annotation_col = as.data.frame(colData(dds))[,1:4],
         annotation_colors = degColors(colData(dds)[1:4],
                                       con_values = c("white",
                                                      "red")
         )
)

data(humanGender)
library(DESeq2)
idx <- c(1:10, 75:85)
dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
                              colData(humanGender)[idx,], design=~group)
degPCA(log2(counts(dse)+0.5), colData(dse),
       condition="group", name="group", shape="group")

#----------------------------------------------------------------
gexprs <- read.delim("fragmentspergene.fragments_per_gene.tsv", 
                     header = T, stringsAsFactors = F)

excl <- c("__alignment_not_unique", "__ambiguous", "__no_feature", 
          "__not_aligned", "__too_low_aQual")
gexprs <- gexprs[!gexprs$feature %in% excl, ]
rownames(gexprs) <- gexprs$feature
gexprs$feature <- NULL

pheno <- read.delim("/FullPheno_RNASeq.txt", header = T)
all(colnames(gexprs) %in% pheno$Samples)
all(colnames(gexprs) == pheno$Names)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=gexprs,
                              colData=pheno, design = ~Patient +Status)

dds <- DESeq(dds)
res <- results(dds)
counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))
degCheckFactors(counts[, 1:6])

degQC(counts, design[["group"]], pvalue = res[["pvalue"]])

th <- HeatmapAnnotation(df = colData(dds),
                        col = degColors(colData(dds), TRUE))

Heatmap(log2(counts(dds) + 0.5)[1:10,],
        top_annotation = th)

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Status"))

vsdT <- vst(dds, blind = T)
vsdF <- vst(dds, blind = F)
#it uses the design formula to calculate the within-group variability (if blind=FALSE) or the across-all-samples variability (if blind=TRUE).
vst_matT <- assay(vsdT)
pcaT <- prcomp(t(vst_matT))

vst_matF <- assay(vsdF)
pcaF <- prcomp(t(vst_matF))


# Create data frame with metadata and PC3 and PC4 values for input to ggplot
library(ggplot2)
df <- cbind(pheno, pcaT$x)
ggplot(df) + geom_point(aes(x=PC1, y=PC2, color = Status))

dfF <- cbind(pheno, pcaF$x)

plotData <- data.frame(pheno, pcaF$x[,1:3])
percentVar <- pcaF$sdev^2/sum(pcaF$sdev^2)

p <- ggplot(plotData, aes(x = PC1, y = PC2, colour = Status))
p + geom_point() + geom_text(aes(label = Samples, hjust = 0, size = 0.5))

plotData <- data.frame(pheno, pcaF$x[,1:3])
p <- ggplot(plotData, aes(x = PC1, y = PC2, colour = Status))
p + geom_point(size = 3) +
  scale_color_manual(values = c("salmon", "navy"))+
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  ggtitle("PCA cartilage non-outilers samples")+
  theme(plot.title = element_text(family = "Helvetica", color="#666666", face="bold", 
                                  size=16, hjust=0.5))

plotData <- data.frame(pheno, pcF$x[,1:3])
p <- ggplot(plotData, aes(x = PC1, y = PC2, colour = Status))
p + geom_point(size = 3) +
  scale_color_manual(values = c("salmon", "navy"))+
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
  geom_text(aes(label = pheno$Samples, hjust = 0.5,vjust = 2, size = 1), cex = 2.5)

#Add lines on the paired samples
p + geom_point(size = 3) +
  scale_color_manual(values = c("salmon", "navy"))+
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
  geom_line(aes(group =pheno$Patient))+
  geom_text(aes(label = pheno$Samples, hjust = 0.5,vjust = 2, size = 1), cex = 2.5)

#Correct for batch effects
mat <- assay(vsdT)
mat <- limma::removeBatchEffect(mat, vsdT$Batch)
assay(vsdT) <- mat
plotPCA(vsdT, intgroup = c("Batch") )
pcaT2 <- prcomp(t(mat))
dfT2 <- cbind(pheno, pcaT2$x)
ggplot(dfT2) + geom_point(aes(x=PC1, y=PC2, color = Batch))

matF <- assay(vsdF)
matF <- limma::removeBatchEffect(matF, vsdF$Batch)

pcaF <- prcomp(t(matF))
dfF2 <- cbind(pheno, pcaF$x)
ggplot(dfF2) + geom_point(aes(x=PC1, y=PC2, color = Batch))

degPCA(vsdT), colData(dds),
       condition="Status", name="Samples", shape="Batch")

p + geom_point() + theme_bw() + geom_text(aes(label = SAMPLES, hjust = 0, size = 1))

## Compute pairwise correlation values
vst_cor <- cor(mat)  
pheatmap(vst_cor)
heat.colors <- brewer.pal(6, "Blues")
pheatmap(vst_cor, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)
