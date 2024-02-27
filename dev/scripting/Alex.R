# Loading the data and packages ---------------------------------

install.packages("Matrix")
install.packages("umap")
install.packages("progress")
install.packages("Seurat")
install.packages(c("Rcpp", "ggplot2", "irlba", "Matrix", "splines", "RColorBrewer", 
                   "dplyr", "bipartite", "pracma", "Rtsne", "R.utils"))
library(Seurat)

library(Matrix)
library(progress)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(umap)


source('dev/utilities/dataprocessingHelperFunctions.R') 
#Load the single cell data for seed 12
data = load('data/GSE226097_seedling_12d_230221.RData')
#loading in the og AraBoxcis network 
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)


# What is contained in the data set?  ---------------------------
print(data)
#the data contains two variables, gbox and clust. 
#gbox is a matrix (rows and columns)
#clust is a vectors (one single column)

#dim give the number of rows and columns 
dim(gbox)

length(clust)

#print names of rows in gbox
rownames(gbox)[1:10]

#print names of cols in gbox
colnames(gbox)[1:10]

#printing part of the matrix gbox
as.matrix(gbox[1:20, 1:3])

#Why are there so many zeros in RNA data? 
#There are two types of zeros:
#  a) Biological zeros: These are actual zeros that signify the 
#     absence of a gene's mRNA.
#  b) Non-biological zeros: These reflect the loss of information due 
#     to inefficiencies throughout the experimental process. 

# Clust contiains the cluster designation of each cell (0-12)
# The clusters were identified using the Seurat Package
clust[1:8]

# We can make a table of the clusters. 
table(clust)

#Making a frequency plot
#Create the barplot
frequencyplot_clustering_seedling_d12 <- plot(table(sort(clust)),
                                        xlab = 'Cluster Name',
                                        ylab = 'Number of Cells',
                                        main = 'Seedling - Day 12')


# Araboxcis -----------------------------------------------------

dim(araboxcis)
araboxcis[1:4,]
hist(araboxcis[,3])

#getting a list of tfs
tfs = unique(araboxcis[,1])
tfSubs = tfs[which(tfs %in% rownames(gbox))]
length(tfSubs) #153

#filtering genes and cells
dim(gbox)
thresh = 0.01
numberGenesPerCell = apply(gbox, 2, function(i){length(which(i > 0))})
includeCells = which(numberGenesPerCell > (thresh*dim(gbox)[1]))
gbox_filtered = gbox[,includeCells]
dim(gbox_filtered)

# Vis using UMAP ------------------------------------------------

gbox.umap <- umap(gbox_filtered, n_neighbors = 10)
png("dev/figures/umap_seedling_d12_k10.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], 
     gbox.umap$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'UMAP Seedling d12, k = 10', 
     xlab = 'UMAP Component 1',
     ylab = 'UMAP Component 2',
     xlim = c(-10, 10),
     ylim = c(-7.5, 7.5))
dev.off()

gbox.umap_15 <- umap(gbox_filtered, n_neighbors = 15)

colours = rainbow(length(unique(clust)))
plot(gbox.umap_15$layout[,1], 
     gbox.umap_15$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'UMAP Seedling d12, k = 15', 
     xlab = 'UMAP Component 1',
     ylab = 'UMAP Component 2')



gbox.umap_30 <- umap(gbox_filtered, n_neighbors = 30)
png("dev/figures/umap_seedling_d12_k30.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.umap_30$layout[,1], 
     gbox.umap_30$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'UMAP Seedling d12, k = 30', 
     xlab = 'UMAP Component 1',
     ylab = 'UMAP Component 2',
     xlim = c(-10, 10),
     ylim = c(-7.5, 7.5))
dev.off()


gbox.umap_100 <- umap(gbox_filtered, n_neighbors = 100)
png("dev/figures/umap_seedling_d12_k100.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.umap_100$layout[,1], 
     gbox.umap_100$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'UMAP Seedling d12, k = 100', 
     xlab = 'UMAP Component 1',
     ylab = 'UMAP Component 2',
     xlim = c(-10, 10),
     ylim = c(-7.5, 7.5))
dev.off()


pca = prcomp(gbox_filtered, scale. = T, rank. = 5)

gbox.pca.umap_10 <- umap(pca$x, n_neighbors = 10)
png("dev/figures/umap_PCA_seedling_d12_k10.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_10$layout[,1],
     gbox.pca.umap_10$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, k = 10', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2',
    xlim = c(-10, 10),
    ylim = c(-7.5, 7.5))
dev.off()

gbox.pca.umap_15 <- umap(pca$x, n_neighbors = 15)
png("dev/figures/umap_PCA_seedling_d12_k15.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_15$layout[,1],
     gbox.pca.umap_15$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, k = 15', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2',
     xlim = c(-10, 10),
     ylim = c(-10, 10))
dev.off()

gbox.pca.umap_30 <- umap(pca$x, n_neighbors = 30)
png("dev/figures/umap_PCA_seedling_d12_k30.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_30$layout[,1],
     gbox.pca.umap_30$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, k = 30', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2',
     xlim = c(-10, 10),
     ylim = c(-10, 10))
dev.off()

gbox.pca.umap_100 <- umap(pca$x, n_neighbors = 100)
png("dev/figures/umap_PCA_seedling_d12_k100.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_100$layout[,1],
     gbox.pca.umap_100$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, k = 100', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2',
     xlim = c(-10, 10),
     ylim = c(-10, 10))
dev.off()



pca = prcomp(gbox_filtered, scale. = T, rank. = 5)

gbox.pca.umap_10 <- umap(pca$x, n_neighbors = 10)
png("dev/figures/umap_PCA_seedling_d12_r5.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_10$layout[,1],
     gbox.pca.umap_10$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, rank = 5', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')
dev.off()

pca = prcomp(gbox_filtered, scale. = T, rank. = 10)

gbox.pca.umap_10 <- umap(pca$x, n_neighbors = 10)
png("dev/figures/umap_PCA_seedling_d12_r10.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_10$layout[,1],
     gbox.pca.umap_10$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, rank = 10', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')
dev.off()

pca = prcomp(gbox_filtered, scale. = T, rank. = 50)

gbox.pca.umap_10 <- umap(pca$x, n_neighbors = 10)
png("dev/figures/umap_PCA_seedling_d12_r50.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_10$layout[,1],
     gbox.pca.umap_10$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, rank = 50', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')
dev.off()

pca = prcomp(gbox_filtered, scale. = T, rank. = 100)

gbox.pca.umap_10 <- umap(pca$x, n_neighbors = 10)
png("dev/figures/umap_PCA_seedling_d12_r100.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_10$layout[,1],
     gbox.pca.umap_10$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, rank = 100', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')
dev.off()

pca = prcomp(gbox_filtered, scale. = T)

gbox.pca.umap_10 <- umap(pca$x, n_neighbors = 10)
png("dev/figures/umap_PCA_seedling_d12_rud.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_10$layout[,1],
     gbox.pca.umap_10$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, rank undefinded', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')
dev.off()


# Attempting Pseudotime ---------------------------------------

