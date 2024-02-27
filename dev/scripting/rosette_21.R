#install.packages('Matrix')
library(Matrix)

#if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
source('dev/utilities/dataprocessingHelperFunctions.R')

#load file
a=load('data/rosette_21d.RData')

#Load the original AraBOXcis network that was trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#check what files are contained in the RNA-seq file
print(a)

#check the dimensions of a table/matrix
dim(gbox)

#check length of a vector
length(clust)

#print first 10 rownamwa of the gbox sequence
rownames(gbox)[1:10]

#print first 10 column names of the gbox
colnames(gbox)[1:10]

#print the first 50 rows and 3 columns
as.matrix(gbox[1:50, 1:3])

#print the cluster designation of the first 8 cells
clust[1:8]

#order the number of cells in each cluster
table(clust)

#plot the table as a graph
plot(table(sort(clust)), xlab='Cluster name', ylab='Number of cells', main='Rosette')

#check dimensions of gboxNetwork22C
dim(araboxcis)

#print the first 4 rows
araboxcis[1:4,]

#create a histogram of the scores
hist(araboxcis[,3])

#list of transcription factors
tfs=unique(araboxcis[,1])

#filter to only have TFs that are in single cell RNA-seq data
tfSubs=tfs[which(tfs %in% rownames(gbox))]

#print the number of TFs
length(tfSubs)

#Unique TFs Rosette: 145

#Check the number of rows and columns of the gbox
dim(gbox)

#Remove cells that have less that 1% of the genes expressed
thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

#Filter to remove genes that are expressed in less than 1% of the cells
numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)


#load package for visualisation
library(umap)

#Let us visualise it using UMAP
gbox.umap <- umap(gbox_filtered)

#Create plot of UMAP
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1],
     gbox.umap$layout[,2],
     col=colours[clust[includeCells]],
     pch=20, main='UMAP Rosette', xlab='UMAP Component 1', ylab='UMAP Component 2')

#PCA
pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], 
     gbox.pca.umap$layout[,2], 
     col=colours[clust[includeCells]], 
     pch=20, main='PCA UMAP Rosette',
     xlab='UMAP Component 1',
     ylab='UMAP Component 2')

#cluster of interest 
cluster_of_interest <- 16

#Visualise whether there are cell type-specific patterns 
colours <- rep('grey', length(unique(clust)))
colours[cluster_of_interest] <- "purple"
plot(gbox.pca.umap$layout[,1],
     gbox.pca.umap$layout[,2],
     col=colours[clust[includeCells]],
     pch=20, main='PCA UMAP Rosette',
     xlab='UMAP Component 1', ylab='UMAP Component 2')


#load package for T-sne
library(Rtsne)
library(ggplot2)

#Create matix and perform T-sne
#Priya helped me with this.
gbox_filtered <- as.matrix(gbox_filtered)
gbox_matrix <- as.matrix(gbox_filtered)

#Specify parameters

tsne_result <- Rtsne(gbox_filtered, prerplexity = 50, theta = 0.5, dims =2)
tsne_result <- Rtsne(gbox_matrix, dims = 2, perplexity = 45)

#Plot T-sne
plot(tsne_result$Y, col=colours[clust[includeCells]], pch = 20, main = "t-SNE plot Rosette_21d")


