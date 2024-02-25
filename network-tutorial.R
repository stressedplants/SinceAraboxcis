library(Matrix)
library(GENIE3)
source('dev/utilities/dataprocessingHelperFunctions.R')

#load data
#Load the single cell data for flowers.
a=load('data/GSE226097_flower_230221.RData')
#Load the original AraBOXcis network that was trained on bulk RNA-seq in seedlings 
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#what is contained in the data set?
#It contains two variables: gbox and clust. gbox is a Matrix and clust is a vector. That means that gbox is a data table with rows and columns, while clust is a single column of data.
print(a)
#check the dimensions of the table/matrix (number of rows and number of column in the table)
dim(gbox)
#check the length of a vector with the function length.
length(clust)
#The rownames of gbox are the names of the genes in Arabidopsis thaliana. print the first 10 rownames.
rownames(gbox)[1:10]
#The columns of gbox are the individual cells.print out the first 10 column names: 
colnames(gbox)[1:10]
#The values in the Matrix are the gene expression values.look at the first 50 rows and the first 3 columns:
as.matrix(gbox[1:50, 1:3])
#clust contains the cluster designation of each cellhe length of clust is the same as the number of columns of gbox,because every cell has a cluster designation. In that paper, they also try to label each cluster with its cell type on the basis of the gene expression pattern of known cell type markers. We can look at that mapping later when we try to interpret the network. print the cluster designations of the first 8 cells:
clust[1:8]
# table to count the number of cells in each cluster:
table(clust)
#plot 
plot(table(sort(clust)), xlab='cluster name', ylab='number of cells')
#check dimensions 
dim(araboxcis)
#look at the first 4 rows.The first column contains the transcription factors. The second column contains the target genes. The third column contains the score of the regulatory edge.select a threshold (often 0.05) for deciding if an edge should be part of the network or not 
araboxcis[1:4,]
#make a histogram of the scores 
hist(araboxcis[,3])
#needed for the next analysis 
tfs=unique(araboxcis[,1])
#filter this to only include transcription factors that are also in the single cell RNA-seq data set. 
tfSubs=tfs[which(tfs %in% rownames(gbox))]
#find out the number of transcription factors 
length(tfSubs)

#First glimpse at Genie3
subGbox=as.matrix(gbox[,which(clust==12)])
net=GENIE3(subGbox, regulators = tfSubs, nTrees=100)
ginieOutput=convertToAdjacency(net, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]