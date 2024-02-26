library(Matrix)

#function source will run all code from a different file
source('dev/utilities/dataprocessingHelperFunctions.R')


#loading with part i'm looking at
a = load('data/GSE226097_rosette_30d_230221.RData')

#Load the original AraBOXcis network that was trained on bulk RNA-seq in seedlings
#set the parameter header as true if you have column names
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)

print(a)

#prints out 2 variables gbox + clust
#gbox is a matrix: a data table with rows and columns
#clust is a vector: single column of data

#function dim tells u no. of rows + columns
dim(gbox)

#function length checks length of cluster
length(clust)

#the rownames of gbox are gene names
#this dataset already been filtered to only have genes that are in the G-box network
#so these genes are either G-box binding TFs or genes with perfect G-box sequence motifs in their promoters

#looks at first 10 rownames in gbox
rownames(gbox)[1:10]

#the output of rownames is a vector
#to access data in a vector, u can specify the id(number) or ids(vector) u want to access within square bracket
# such as vec[id] or vec[ids]

#the column names of gbox are the individual cells; they all have unique ids

colnames(gbox)[1:10]

#the values in the matrix gbox are gene expression values
#there r many 0s bcos 
as.matrix(gbox[1:50, 1:3])
#looks at first 50 rows, and first 3 columns
#as.matrix turns a Matrix (with capital M) into a normal matrix, which is a normal table in R
#to get submatrix of Matrix u do mat[rows, cols]
#if u want certain rows + all columns  mat[rows, ]
#if u want all rows + certain columns mat[,cols]

#the variable clust contains cluster designation of each cell
#the clusters were identified using Seurat package

#prints cluster designations of first 8 cells
clust[1:8]

#function table counts number of cells in each cluster
table(clust)

#plot this as a graph
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Rosette 30 Days')

#---------------------------------------------------------------------
#looking at gboxNetwork22C.csv
dim(araboxcis)

#looks at first 4 rows
araboxcis[1:4,]

#column 1 is the TF
#column 2 is the target gene
#column 3 is the score of the regulatory edge

#this networks includes all edges with score>0, but usually u would choose a threshold
#normally 0.05 to decide if u want to include an edge or not in the network

#makes a histogram of the scores, Scores is column 3 + u want all rows
hist(araboxcis[,3])

#gets a list of TFs
tfs=unique(araboxcis[,1])
#unique function creates a new vector with all duplicates removed

#to filter this to only include TFs that are also in the scRNA-seq data
#some TFs are expressed so lowly that they aren't observed in scRNA-seq

tfSubs=tfs[which(tfs %in% rownames(gbox))]

#find out how many TFs there are
length(tfSubs)

#---------------------------------------------------------------------
#Filter Genes and cells with very low values

#gives rows + columns of gbox
dim(gbox)

#get rid of cells with less than 1% of the genes expressed at all

thresh=0.01

#calculates no. of genes for each cell in gbox
#function counts no. of genes in each cell where count is bigger than 0
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})

#will only include the cells which meet the threshold for the no. of genes per cell
#so it multiplies the threshold by the total no. of gene (dim(gbox)[1])
#dim(Gbox)[1] gives you the no. of rows which = no. of genes
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

#makes new matrix called gbox_Filtered which filters cells in gbox where the
#gene expression meets threshold
gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

#get rid of genes expressed in less than 1% of the cells

numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))

gbox_filtered=gbox_filtered[includeGenes,]

dim(gbox_filtered)

#Gbox filtered when thresh=0.01 = 1520 37816
#Gbox filtered when thresh=0.05 = 784 13865 - a lot of points taken away but no clearer clusters

#---------------------------------------------------------------------
#Visualisation with UMAP

#high dimensionl data -> 2D data with UMAP

install.packages('umap')
library(umap)

#UMAP IS VERY SLOW - can take like 10m
gbox.umap <- umap(gbox_filtered)

#do the cell type clusters group together if we only look at G-box related genes?
#pch=20 means filled circles; used to specify plotting symbol
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Rosette 30d', xlab='UMAP Component 1', ylab='UMAP Component 2')

#---------------------------------------------------------------------
#PCA before UMAP

pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Rosette 30 Days', xlab='UMAP Component 1', ylab='UMAP Component 2')

#---------------------------------------------------------------------
#trying to. just colour one cluster
cluster_of_interest <- 1

#colours vectors is set to grey + the cluster of interest is coloured red
#rep: repeat for all clusters

colours <- rep('grey', length(unique(clust)))
colours[cluster_of_interest] <- 'red'
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Rosette 30 Days - Cluster 1 in red', xlab='UMAP Component 1', ylab='UMAP Component 2')

#---------------------------------------------------------------------
#trying to do t-SNE

install.packages("Rtsne")
library(Rtsne)


gbox_filtered <- as.matrix(gbox_filtered)

tsne_result <- Rtsne(gbox_filtered)

tsne_result <- Rtsne(gbox_filtered, perplexity = 50, theta = 0.5, dims = 2)

#the Y are the matrix of t-SNE coordinates made by the Rtsne function
colours=rainbow(length(unique(clust)))
plot(tsne_result$Y, col=colours[clust[includeCells]], pch = 20, main = "t-SNE Plot")


#---------------------------------------------------------------------
