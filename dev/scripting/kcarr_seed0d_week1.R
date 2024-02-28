 
# load all packages:
library(monocle3)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(umap)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(Matrix)
library(patchwork)
library(R.utils)
library(shiny)
library(SingleCellExperiment)
library(mclust)

# load in unique dataset (seed 0d)
a = load('data/GSE226097_seed_0d_230221.RData')

#load araboxcis network: 

araboxcis = read.csv('data/gboxNetwork22C.csv')

print(a)
dim(gbox)
length(clust)

#print the first 10 rownames of gbox

rownames(gbox)[1:10]

#print the first 10 column names

colnames(gbox)[1:10]

#look at the first 50 rows and 3 columns to examine the 0's

as.matrix(gbox[1:50, 1:3])

#print the cluster designations of the first 8 cells in clust

clust[1:8]

#use the table function to count the number of cells in each 
# cluster

table(clust)

#sort the cluster info into a table 

plot(table(sort(clust)), xlab = 'cluster name', ylab = 'number of cells', main = 'seed 0d')

dim(araboxcis)

araboxcis[1:4,]

#the first column contains transcription factors, second column is target genes, 
# third column is the score of the regulatory edge. 

hist(araboxcis[,3])

#need a list of transcription factors for next week

tfs=unique(araboxcis[,1])

#filter tfs to include only tfs that are also in the ssRNA-Seq set. 

tfSubs = tfs[which(tfs %in% rownames(gbox))]

length(tfSubs)

#there are 153 unique tf's that are in both datasets

#Seed0d: 153

#Filter genes and cells with very low values: 

#get rid of cells that have less than 1% of genes expressed at all 

thresh = 0.01
numberGenesPerCell = apply(gbox, 2, function(i){length(which(i>0))})
includeCells = which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

numberCellPerGene = apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes = which(numberCellPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered = gbox_filtered[includeGenes,]
dim(gbox_filtered)

#visualization with UMAP 


library(umap)

gbox.umap <- umap(gbox_filtered)

colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Seed0d', xlab='UMAP Component 1', ylab='UMAP Component 2')

# try the UMAP with diff number of neighbors to see if it changes the appearance

num_neighbors <- 10

gbox_umap_10 <- umap(gbox_filtered, n_neighbors = num_neighbors)

plot(gbox_umap_10$layout[,1], gbox_umap_10$layout[,2], col = colours[clust[includeCells]], pch = 20, main = 'UMAP Seed0d 200', xlab = "UMAP Component 1", ylab = "UMAP Component 2")

#It looks the same with n_neighbors set to 10. 


#try PCA Visualisation too: 

pca = prcomp(gbox_filtered, scale.=T, rank.=200)
gbox.pca.umap <- umap(pca$x, n_neighbors = num_neighbors)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Seed0d', xlab='UMAP Component 1', ylab='UMAP Component 2')

#perform tSNE for further comparison: 

library(Rtsne)
#first, use as.matrix to save gbox_filtered as a matrix for tsne. 
gbox_filtered_matrix <- as.matrix(gbox_filtered)

#next, run the tsne and save it as a variable called tsne_result:
tsne_result <- Rtsne(gbox_filtered_matrix, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

#next, plot the tsne_result to view the results: 
plot(tsne_result$Y, col = "blue", pch = 16, main = "Seed0d tSNE Visualisation")

#the resulting plot is a little difficult to interpret, change the colours and visualisations: 

plot(tsne_result$Y, col = colours[clust[includeCells]], pch = 20, main = "Seed_0d t-SNE Visualisation", xlab = "t-SNE 1", ylab = "t-SNE 2")





#Try using Seurat to perform clustering/UMAP

library(Seurat)



if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages('Seurat')
remotes::install_github('satijalab/seurat-wrappers')
library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)



# Run UMAP in Seurat: 

#nFeature_RNA is the number of genes detected in each cell, nCount_RNA is the number of total molecules detected in a cell. 
# this is a QC step. 
Seurat_Object_UMAP <- CreateSeuratObject(counts = gbox_filtered, project = "Araboxcis UMAP")
meta.data <- Seurat_Object_UMAP@meta.data
head(meta.data)

#this graph is 
FeatureScatter(Seurat_Object_UMAP, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', group.by = 'orig.ident')

# Normalize the data
Seurat_Object_UMAP <- NormalizeData(Seurat_Object_UMAP, normalization.method = "LogNormalize", scale.factor = 10000)
Seurat_Object_UMAP_simple <- CreateSeuratObject(counts = gbox_filtered, project = "Araboxcis UMAP2")
Seurat_Object_UMAP_simple <- NormalizeData(Seurat_Object_UMAP_simple)

#identify highly variable features: 
Seurat_Object_UMAP <- FindVariableFeatures(Seurat_Object_UMAP, selection.method = "vst", nfeatures = 2000)
top10_attemp2 <- head(VariableFeatures(Seurat_Object_UMAP), 10)

#Plot variable features
plot1_attempt2 <- VariableFeaturePlot(Seurat_Object_UMAP)
LabelPoints(plot = plot1_attempt2, points = top10_attemp2, repel = T)

#Scaling Data
all.genes2 <- rownames(Seurat_Object_UMAP)
Seurat_Object_UMAP <- ScaleData(Seurat_Object_UMAP, features = all.genes2)

#Perform Dimensionality Reduction
Seurat_Object_UMAP <- RunPCA(Seurat_Object_UMAP, features = VariableFeatures(object = Seurat_Object_UMAP))

#Determine Dimensionality

ElbowPlot(Seurat_Object_UMAP)

#Clustering data
Seurat_Object_UMAP <- FindNeighbors(Seurat_Object_UMAP, dims = 1:20)
Seurat_Object_UMAP <- FindClusters(Seurat_Object_UMAP, resolution = 0.5)
# found 9 clusters

DimPlot(Seurat_Object_UMAP, group.by = "RNA_snn_res.0.5", label = T)
# Set identity of clusters
Idents(Seurat_Object_UMAP) <- "RNA_snn_res.0.5"

#Perform UMAP!
Seurat_Object_UMAP <- RunUMAP(Seurat_Object_UMAP, dims = 1:20)

DimPlot(Seurat_Object_UMAP, reduction = 'umap')

#try tSNE: 

Seurat_Object_UMAP <- RunTSNE(object = Seurat_Object_UMAP, dims = 1:10, do.fast = T)
DimPlot(Seurat_Object_UMAP, reduction = 'tsne')

#I found that the Seurat method seemed to identify clusters, whereas the 
# non-Seurat method seemed to show a pseudotime trajectory. I'm not sure 
# why this is. 