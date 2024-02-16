#Loading data
# - Installing the matrix package 
library(Matrix)
# - Adding source 
source("C:/Users/tammy/OneDrive/Documents/GitHub/SinceAraboxcis/dev/utilities/dataprocessingHelperFunctions.R")
# - Loading the dataset (seedling day 6)
a = load("/Users/tammy/OneDrive/Documents/GitHub/SinceAraboxcis/data/GSE226097_seedling_6d_230221.RData")
# - Loading Araboxcis network 
araboxcis = read.csv("/Users/tammy/OneDrive/Documents/GitHub/SinceAraboxcis/data/gboxNetwork22C.csv", header=T)

#Inspecting the dataset
# - Seeing what variables were saved
print(a)
# - Checking dimensions of gbox matrix 
dim(gbox)
# - Checking length of vector 
length(clust)
# - Displaying first 10 row names
rownames(gbox)[1:10]
# - Displaying first 10 columns
colnames(gbox)[1:10]
# - Displaying first 50 rows and first 3 columns 
as.matrix(gbox[1:50, 1:3])
# - Displaying the cluster designations for first 8 cells 
clust[1:8]
# - Identifying number of cells per cluster 
table(clust)

#Plotting clusters 
plot(table(sort(clust)),
     xlab = "cluster name", 
     ylab = "number of cells", 
     main = "Seedling6d")

#Inspecting gboxNetwork22C datafile 
# - Checking dimensions of matrix 
dim(araboxcis)
# - Displaying first 4 rows 
araboxcis[1:4,]

#Plotting histogram of araboxcis
# - Histogram of 3rd column (score of regulatory edge)
hist(araboxcis[,3])

#Transcription factors 
# - List transcription factors as variable
tfs = unique(araboxcis[,1])
# - Find only transcription factors that are in the single cell RNA-seq dataset
tfSubs = tfs[which(tfs %in% rownames(gbox))]
# - Length of transcription factors 
length(tfSubs)

#Filtering genes and cells with very low values 
# - Removing cells with <1% genes
# - Defining threshold value 
thresh = 0.01  
# - Defining number of genes per cell
numberGenesPerCell = apply(gbox, 2, 
                           function(i){length(which(i>0))})
# - Choosing which cells to filter out 
includeCells = which(numberGenesPerCell>(thresh*dim(gbox)[1]))
# - Filtering the gbox 
gbox_filtered = gbox[,includeCells]
# - Checking the new dimensions after filtering 
dim(gbox_filtered)
# - Removing genes expressed in <1% of cells 
# - Defining number of cells per gene 
numberCellsPerGene = apply(gbox_filtered, 1, 
                           function(i){length(which(i>0))})
# - Choosing which genes to filter out 
includeGenes = which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
# - Filtering the gbox 
gbox_filtered = gbox_filtered[includeGenes,]
# - Checking the new dimensions after filtering 
dim(gbox_filtered)

#Visualisation with UMAP 
# - Loading UMAP library 
install.packages('umap')
library(umap)
# - Visualising gbox using UMAP 
gbox.umap <- umap(gbox_filtered)
# - Setting colours
colours = rainbow(length(unique(clust)))
# - Making the plot 
plot(gbox.umap$layout[,1], 
     gbox.umap$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'UMAP Seedlingd6', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')

#Visualisation with UMAP & PCA
# - Defining PCA 
pca = prcomp(gbox_filtered, scale.=T, rank.=5)
# - Visualising gbox using PCA
gbox.pca.umap <- umap(pca$x)
# - Setting colours 
colours = rainbow(length(unique(clust)))
# - Making the plot 
plot(gbox.pca.umap$layout[,1], 
     gbox.pca.umap$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'PCA UMAP Seedlingd6', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2') 