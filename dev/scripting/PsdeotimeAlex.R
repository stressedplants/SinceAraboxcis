library(GENIE3)

library(Seurat)

library(Matrix)
library(progress)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(umap)

source("http://bioconductor.org/biocLite.R")
install("qlcMatrix")
library(monocle)

source('dev/utilities/dataprocessingHelperFunctions.R') 
#Load the single cell data for seed 12
data = load('data/GSE226097_seedling_12d_230221.RData')
#loading in the og AraBoxcis network 
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)



# Araboxcis -----------------------------------------------------

#getting a list of tfs
tfs = unique(araboxcis[,1])
tfSubs = tfs[which(tfs %in% rownames(gbox))]
length(tfSubs) #153

#filtering genes and cells
dim(gbox)
thresh = 0.06 #change to one later
numberGenesPerCell = apply(gbox, 2, function(i){length(which(i > 0))})
includeCells = which(numberGenesPerCell > (thresh*dim(gbox)[1]))
gbox_filtered = gbox[,includeCells]
dim(gbox_filtered)

# Vis using UMAP ------------------------------------------------

gbox.umap <- umap(gbox_filtered)
#png("dev/figures/umap_seedling_d12_k10.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], 
     gbox.umap$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'UMAP Seedling d12', 
     xlab = 'UMAP Component 1',
     ylab = 'UMAP Component 2',
     xlim = c(-10, 10),
     ylim = c(-7.5, 7.5))
#dev.off()

pca = prcomp(gbox_filtered, scale. = T, rank. = 5)

gbox.pca.umap_10 <- umap(pca$x, n_neighbors = 10)
#png("dev/figures/umap_PCA_seedling_d12_k10.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_10[,1],
     gbox.pca.umap_10[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, k = 10', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')
  
#dev.off()

#attempt to cluster
pca <- prcomp(gbox_filtered, scale. = TRUE, rank. = 5)

# Obtain UMAP embeddings
library(uwot)
umap_result <- umap(pca$x, n_neighbors = 10)

# Cluster cells using K-means clustering
num_clusters <- 10  # Number of clusters you want to generate
cluster_labels <- kmeans(umap_result, centers = num_clusters)$cluster

# make a suerate object ---------------------------------------------
library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(patchwork)
library(dplyr)
library(Signac)
library(Matrix)
remove.packages("irlba")
install.packages("irlba", type = "source", force = TRUE)
library(irlba)


source('dev/utilities/dataprocessingHelperFunctions.R') 
#Load the single cell data for seed 12
data = load('data/GSE226097_seedling_12d_230221.RData')
#loading in the og AraBoxcis network 
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)

tfs = unique(araboxcis[,1])
tfSubs = tfs[which(tfs %in% rownames(gbox))]
length(tfSubs) #153

#filtering genes and cells
dim(gbox)
thresh = 0.01 #change to one later
numberGenesPerCell = apply(gbox, 2, function(i){length(which(i > 0))})
includeCells = which(numberGenesPerCell > (thresh*dim(gbox)[1]))
gbox_filtered = gbox[,includeCells]
dim(gbox_filtered)


gbox_matrix_filtered <- as.matrix(gbox_filtered)
cds <- new_cell_data_set(gbox_matrix_filtered)
cds_preprocessed <- preprocess_cds(cds, num_dim = 5)
plot_pc_variance_explained(cds_preprocessed)
cds_reduce_dim <- reduce_dimension(cds_preprocessed)

cds_cluster <- cluster_cells(cds_reduce_dim, reduction_method = ("UMAP"), k = 15)
cds_cluster_learnt <- learn_graph(cds_cluster)
plot_cells(cds_cluster_learnt)
plot_cells(cds_cluster_learnt, color_cells_by = "partition")
cds_ordered <- order_cells(cds_cluster_learnt)


plot_cells(cds_ordered,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5)




