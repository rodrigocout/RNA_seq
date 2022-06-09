deg <- function(gexp, pheno, design, condition, dirname=".")
{
  require("DESeq2")
   dds <- DESeqDataSetFromMatrix(countData=gexp,
                         colData=pheno, design = design)
   deg <-  DESeq(dds)
   DESeq_results <- results(deg, contrast=condition, alpha=0.05)
   DESeq_results <- as.data.frame(DESeq_results)
   DESeq_results

write.table(DESeq_results, paste(dirname, "DESeqResults.txt"), 
            sep = "\t", quote = F)
}

DEG <- deg(plasm, info, design = ~Condition, 
           condition = c("Condition","Case", "Control"), 
           dirname = "/RNASeq_pipelines/")

#----------
#VST norm
VstNorm <- function(gexp, pheno, condition = ".", dirname=".", plotpdf = FALSE)
{
require("DESeq2")
dds <- DESeqDataSetFromMatrix(countData=gexp,
                                colData=pheno, design = ~1)
vst <- varianceStabilizingTransformation(dds, blind = T)
MatVST<- assay(vst)
MatVST

pc <- prcomp(t(MatVST))
percenVar <- pc$sdev^2/sum(pc$sdev^2)
plotData <- data.frame(pheno, pc$x[,1:2])
p <- ggplot(plotData, aes(x = PC1, y = PC2, colour = condition))

if(plotpdf)pdf(file=paste(dirname,"PCAPlot.pdf"))
p + geom_point(size = 3) +
  scale_color_manual(values = c("salmon", "blue"))+
  xlab(paste0("PC1: ", round(percenVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percenVar[2] * 100), "% variance")) 
if(plotpdf)dev.off()
}
