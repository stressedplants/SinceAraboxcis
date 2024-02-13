#Ellie seedling3d

#install.packages('Matrix')
library(Matrix)
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
source('dev/utilities/dataprocessingHelperFunctions.R')

#Load the single cell data for seedling3d.
a=load('data/GSE226097_seedling_3d_230221.RData')

#Load the original AraBOXcis network that was trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#what is contained in the data set?
print (a)
#gbox is a matrix (data table with rows and columns) and clust is a vector (single column of data). 
#check the dimensions of matrix and the length of clust
dim(gbox)
length(clust)
#The rownames of gbox are the names of the genes in Arabidopsis thaliana. This 
#dataset has already been filtered to only contain the genes that are in the G-box network. 
#This means that these genes either are G-box binding transcription factors or genes
#with perfect G-box sequence motifs in their promoters
rownames(gbox)[1:10]
colnames(gbox)[1:10]
#will print the 10 names of the rows and columns 

#values in the matrix are gene expression values
#view the first 50 rows and 3 columns
as.matrix(gbox[1:50, 1:3])

#clust contains the cluster designation of each cell. the clusters were identified using seurat
#print the cluster designations of the first 8 cells 
clust[1:8]
#use the function table to count the number of cells in each cluster 
table(clust)
#this can be plotted as a graph, using plot, bars can be sorted by height using the sort function
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Seedling3d')

#check the dimensions of the gboxNetwork22C (also matrix)
dim(araboxcis)
#look at the first 4 rows 
araboxcis[1:4,]
#1st col= transcription factors, 2nd col= target genes. 3rd col= score of the regulatory edge.
#make a histogram of the scores
hist(araboxcis[,3])

#list of transcription factors with a new vector and duplicates removed 
tfs=unique(araboxcis[,1])
#filter to only include transcription factors that are also in the single cell RNA-seq data set.
tfSubs=tfs[which(tfs %in% rownames(gbox))]
#how many transcription factors?
length(tfSubs)

#Filter genes and cells with very low values 
#get rid of cells that have less than 1% of the genes expressed at all.
dim(gbox)
thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)
#get rid of genes that are expressed in less than 1% of the cells
numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)

First visualisation with UMAP
#install.packages('umap')
library(umap)
#visualise it using UMAP
gbox.umap <- umap(gbox_filtered)
#Do the cell type clusters group together if we only look at G-box related genes?
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Seedling_3d', xlab='UMAP Component 1', ylab='UMAP Component 2')

#PCA UMAP 
pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Seedling_3d', xlab='UMAP Component 1', ylab='UMAP Component 2')