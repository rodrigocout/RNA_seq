suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(ggplot2)
  library(scran)
  library(batchelor)
  library(BiocSingular)
  library(BiocNeighbors)
  library(Rtsne)
  library(rsvd)
  library(tidyr) 
  library(dplyr)
  library(igraph)
})

## Download the data and set row names to gene symbols whenever possible
sce <- readRDS(gzcon(url("https://ivanek.github.io/analysisOfGenomicsDataWithR/data/10/SCE_MammaryEpithelial_x3.rds?raw=true")))
rownames(sce) <- scater::uniquifyFeatureNames(
  ID = rownames(sce), 
  names = as.character(rowData(sce)$gene.symbols)
)
#Subsample cells to speed up processing:
set.seed(42)
n=3000
sce <- sce[, sample(1:ncol(sce), n )  ]
## Dataset compostion per cell type and study:  
table(colData(sce)$study , colData(sce)$cell.class)


## Check library size amd  number of detected gene distribution in the three datasets:
ggplot(data.frame(colData(sce)), 
       aes(x = study, y = library_size, color=study)) + 
  geom_violin() + theme_bw() + 
  ylab("Total UMI counts per cell") + 
  ggtitle("Library size distribution per study")

ggplot(data.frame(colData(sce)), 
       aes(x = study, y = detected_genes, color=study)) + 
  geom_violin() + theme_bw() + 
  ylab("NODG") + 
  ggtitle("Number of detected genes per study")


# We first normalize all cells for library size.
assays(sce )[["lognorm"]] <- log2(sweep( counts(sce),2,sce$library_size ,FUN="/")*1e3 +1)
reducedDim(sce, "PCA" )  <- rsvd::rpca(t( assay(sce,"lognorm") ),k=32,retx=TRUE,center=TRUE,scale=FALSE)$x
reducedDim(sce, "TSNE" ) <- Rtsne( reducedDim(sce,"PCA"), perplexity = 30, initial_dims=32, pca=FALSE, theta=0.3)$Y #~5"-20" run time

cowplot::plot_grid(scater::plotTSNE(sce, colour_by = "study" ),
                   scater::plotTSNE(sce, colour_by = "cell.class"))



#Batch effect correction

sce.linCor  <- batchelor::rescaleBatches(sce, batch=sce$study)
colData(sce.linCor) <-  colData( sce )  
reducedDim(sce.linCor, "PCA" )  <- rsvd::rpca( t( assay(sce.linCor,"corrected") ) ,k=32,retx=TRUE,center=TRUE,scale=FALSE)$x
reducedDim(sce.linCor, "TSNE" ) <- Rtsne(reducedDim(sce.linCor, "PCA"), perplexity = 30, initial_dims=32, pca=FALSE, theta=0.3)$Y #~15-60 seconds run time
cowplot::plot_grid(scater::plotTSNE(sce.linCor, colour_by = "study" ),
                   scater::plotTSNE(sce.linCor, colour_by = "cell.class")
)


d <- 32
FMNN.out <-  batchelor::fastMNN( assays(sce)[["lognorm"]],  batch=sce$study , use.dimred="PCA", d=d ) 
# Notice that the object returned is a single cell experiment.
reducedDim (sce, "PCA.FMNN" ) <- reducedDim(FMNN.out, "corrected") 
reducedDim( sce, "TSNE.FMNN" ) <- Rtsne( reducedDim(sce, "PCA.FMNN"), perplexity = 30, initial_dims=64, pca=FALSE,num_threads=32,theta=0.25)$Y
cowplot::plot_grid(scater::plotReducedDim(sce, colour_by = "study", dimred="TSNE.FMNN" ),
                   scater::plotReducedDim(sce, colour_by = "cell.class", dimred="TSNE.FMNN")
)

# We can also check the proportion of variance lost from each batch
# and in each merge step of FMNN. In the following matrix rows
# correspond to merge steps and columns to the three batches (according to input order):
metadata(FMNN.out)$merge.info$lost.var

#Clustering

library(scran)
# Create a graph of shared nearest neighbors: Two cells are connected iff 
g <- scran::buildSNNGraph(sce, k=25, use.dimred = 'PCA') # Build SNN graph
clust <- igraph::cluster_louvain(g)$membership # use the louvain method for community detection
table(clust)

sce$clust <- factor(clust) 
cowplot::plot_grid(scater::plotTSNE(sce, colour_by = "study"),
                   scater::plotTSNE(sce, colour_by = "cell.class"),
                   scater::plotTSNE(sce, colour_by = "clust")
)


#Differential expression and identification of marker genes.
sce.vis <- sce[,sce$study=="vis"] # Extract the cells from one specific study

markers_all <- scran::findMarkers(
  sce.vis, groups = sce.vis$cell.class, 
  test.type = "wilcox", assay.type="lognorm",
  pval.type = "all", direction = "up"
)

head(markers_all[[1]])

# Violi plots of expression for a set of features on defined groups
scater::plotExpression(sce.vis, features = c("Csn3", "Trf"), 
                       x = "cell.class")

#Colour a projection according to an identified marker gene:
cowplot::plot_grid(scater::plotTSNE(sce.vis, colour_by = "Csn3"),
                   scater::plotTSNE(sce.vis, colour_by = "cell.class"))

#Heatmap for the top two marker genes per cell group:
scater::plotHeatmap(sce.vis, features = unique(unlist(lapply(markers_all, function(w) rownames(w)[1:2]))),
                    columns = colnames(sce.vis)[order(sce.vis$cell.class)],
                    colour_columns_by = "cell.class", cluster_cols = FALSE,
                    show_colnames = FALSE, cluster_rows = FALSE)

