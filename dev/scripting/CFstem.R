# Calem Stem
library(Matrix)
library(umap)

# BiocManager already installed

# Alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
source('dev/utilities/dataprocessingHelperFunctions.R')

# Load the single cell data for stem.
# The function load(‘filename.RData’) loads R data objects.
a=load('data/GSE226097_stem_230221.RData')

# a contains two variables: gbox and clust. 
# gbox is a Matrix (values in the Matrix are the gene expression values). 
# gbox is a data table with rows and columns
# rownames of gbox are the names of the genes in Arabidopsis thaliana
# gbox contains only the genes that are in the G-box network
# these genes either are G-box binding transcription factors or genes with perfect G-box sequence motifs in their promoters.
# columns of gbox are the individual cells
# clust is a vector.
# clust is a single column of data
# Each cell in the dataset is assigned to a specific cluster based on its gene expression profile


# Load the original AraBOXcis network that was trained on bulk RNA-seq in seedlings
# read.csv(‘filename.csv’) reads in comma separated tables. parameter header=TRUE if there is column names
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

# Use dim function to the number of rows and number of columns in the table
dim(gbox)

# use length function to check the number of elements of the clust vector.
length(clust)

# print the first 10 rownames of gbox. (genes)
rownames(gbox)[1:10]

# print out the first 10 column names of gbox (cells)
colnames(gbox)[1:10]

# look at the first 50 rows and the first 3 columns
as.matrix(gbox[1:50, 1:3])


# print the cluster designations of the first 8 cells
clust[1:8]

# use table function to count the number of cells in each cluster (it gives a frequency table)
table(clust)

# Open a PNG graphics device
png("/Applications/_Level 7 Bioinformatics/Semester 2/Group Project 69/Figures/stem_cluster_plot.png")

# plot this frequency table as a graph
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Stem')

# Close the graphics device
dev.off()


# data also read from a file named gboxNetwork22C which also contains a matrix.
# check dimensions 
dim(araboxcis)

# view first 4 rows
araboxcis[1:4,]

# 1st col = transcription factors
# 2nd col = target genes
# 3rd col = score of the regulatory edge 

# this network is including all the edges that have scores bigger that 0, but usually you would select a threshold (often 0.05) for deciding if an edge should be part of the network or not.
# make a histogram of the scores
hist(araboxcis[,3])

# list the transcription factors
# The unique function creates a new vector with all duplicates removed.
tfs=unique(araboxcis[,1])

# filter to only include transcription factors that are also in the single cell RNA-seq data set
tfSubs=tfs[which(tfs %in% rownames(gbox))]
# ‘vec1 %in% vec2’ produces a vector of length(vec1). The vector contains TRUE wherever an element of vec1 is in vec2, and FALSE otherwise.

# use length function to find out the number of transcription factors.
length(tfSubs) # = stem 139




# filter genes and cells with low values 

# get rid of cells that have less than 1% of the genes expressed at all
thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered) # =2096 17114


# get rid of genes that are expressed in less than 1% of the cells
numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered) # =1645 17114




# PCA plot 

# Open a PNG graphics device
png("/Applications/_Level 7 Bioinformatics/Semester 2/Group Project 69/Figures/PCA_UMAP_stem_plot.png")

pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Stem', xlab='UMAP Component 1', ylab='UMAP Component 2')

# Close the graphics device
dev.off()




# first visualisation with UMAP
# project this high-dimensional data into 2 dimensions using UMAP
gbox.umap <- umap(gbox_filtered)

#Do the cell type clusters group together if we only look at G-box related genes?

# Open a PNG graphics device
png("/Applications/_Level 7 Bioinformatics/Semester 2/Group Project 69/Figures/UMAP_stem_plot.png")

colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Stem', xlab='UMAP Component 1', ylab='UMAP Component 2')

# Close the graphics device
dev.off()


