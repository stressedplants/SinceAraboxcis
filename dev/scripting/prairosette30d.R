library(Matrix)

#function source will run all code from a different file
source('../utilities/dataprocessingHelperFunctions.R')


#loading with part i'm looking at
a = load('../../data/GSE226097_rosette_30d_230221.RData')

#Load the original AraBOXcis network that was trained on bulk RNA-seq in seedlings
#set the parameter header as true if you have column names
araboxcis = read.csv('../../data/gboxNetwork22C.csv', header = T)

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
#seurat

source('../utilities/SeuratUMAP.r')
install.packages("Seurat")
SeuratUMAP_function(gbox_filtered)

#Week 3---------------------------------------------------------------
#finding out the average gene expression of each gene across each cluster

#make clust not be factors bcos rn clust is assumed to be a factor
#turns clust into a numeric vector + assigns it to new variable 'clustAsNumbers'
clustAsNumbers=as.numeric(paste(clust))

#calculates mean expression of the gene within each cluster
geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust))-1), function(j){
    ids=which(clustAsNumbers==j)
    mean(i[ids])
  })
})

#sets column names of the geneExpByCLuster to be the same as the row names
#of the gbox matrix
colnames(geneExpByCluster)=rownames(gbox)

#shows no. of rows + columns of the geneExpByCluster matrix
dim(geneExpByCluster)
#17 2098

#--------------------------------------------------------------------
#Visualising large tables as heatmaps
install.packages('pheatmap')
library('pheatmap')

#each value in the column in geneExpByCLuster wil be scaled 
pheatmap(geneExpByCluster, scale='column')

clustLabs=read.table('../../data/clusterLabels.txt', header=T, sep='\t')

#gets all unique names + assigns samples
unique(clustLabs[,'Organ'])

organ='Rosette 30d'
simpleNames=clustLabs[which(clustLabs[,'Organ']==organ), "Cell.type.suggested"]
print(simpleNames)

rownames(geneExpByCluster)=simpleNames

pheatmap(geneExpByCluster, scale='column', main='Rosette 30d')

write.table(geneExpByCluster, '../../data/Rosette30d_avgExpressionByCluster.txt')

#---------------------------------------------------------------------
#Looking at TFs that bind to perfect Gboxes - to see if these TFs are
#in certain cell types

#gets all unique TFs in araboxcis - TFs in 1st column of araboxcis
tfs=unique(araboxcis[,1])

#gets TFs in araboxcis that is also in gbox
tfSubs=tfs[which(tfs %in% rownames(gbox))]

pheatmap(geneExpByCluster[,tfSubs], scale='column', main='Rosette 30d')

#---------------------------------------------------------------------
#investigating correlation between expression in TFs + potential downstream targets

#to find the correlation between every TF and each potential downstream target across
#the different cell types

#cor calculates coefficient via spearman method
#takes each tf in tfSubs and calculates the correlation between the tf
#and each gene in geneExpByCluster
corMat=sapply(tfSubs, function(tf){
  apply(geneExpByCluster, 2, function(gene){
    cor(geneExpByCluster[,tf], gene, method='spearman')
  })
})

dim(corMat) #2098 144

pheatmap(corMat, main='Rosette 30d')

#---------------------------------------------------------------------
#Looking at correlation between TF expression + downstream target genes
#in individual cells at cell type level

#trying to find TFs that are highly correlated with downstream genes 
# and also if these genes are highly correlated with one another at the
#individual cell level

#filters matrix for where value is >0.8 & is 1 exactly bcos that's 
#just a self correlation - makes it TRUE
id=which(corMat>0.8 & corMat!=1, arr.ind = TRUE)
dim(id) #979 2

dim(corMat)
dim(id)

tVal=apply(id, 1, function(i){
  row=rownames(corMat)[i[1]]
  col=colnames(corMat)[i[2]]
  
  inBoth=length(which(gbox[row,]>0 & gbox[col,]>0))
  inNone=length(which(gbox[row,]==0 & gbox[col,]==0))
  inFirst=length(which(gbox[row,]>0 & gbox[col,]==0))
  inSecond=length(which(gbox[row,]==0 & gbox[col,]>0))
  matTemp=matrix(c(inBoth, inFirst, inSecond, inNone), ncol=2)
  c(inBoth, inFirst, inSecond, inNone, (inBoth*inNone)/(inFirst*inSecond))
})

hist(log(tVal[5,]), main='Rosette 30d', xlab='log odds ratio')

#---------------------------------------------------------------------
#Compare between datasets

thresh=exp(1)
idDoublePositive=which(tVal[5,]>thresh)
doublePositive=id[idDoublePositive,]
doublePositive[,1]=rownames(corMat)[doublePositive[,1]]
doublePositive[,2]=colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive=cbind(doublePositive[,2], doublePositive[,1])
doublePositive=cbind(doublePositive, tVal[5,idDoublePositive])

#save file
write.table(doublePositive, file='../../data/Rosette30d_doublePositives.txt', sep='\t', row.names=F)

thresh=1
idSimpson=which(tVal[5,]<thresh)
Simpson=id[idSimpson,]

Simpson[,1]=rownames(corMat)[Simpson[,1]]
Simpson[,2]=colnames(corMat)[as.numeric(Simpson[,2])]
Simpson=cbind(Simpson[,2], Simpson[,1]) #so TF comes before target
Simpson=cbind(Simpson, tVal[idSimpson])

#save file
write.table(Simpson, file='../../data/Rosette30d_SimpsonPairs.txt', sep='\t', row.names=F)

#---------------------------------------------------------------------
#Making a network from scratch

library(GENIE3)

net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=50)

save(net, file='../../rosette30d_network_nTree_50.RData')#check tree number

#convert file to a table

ginieOutput=convertToAdjacency(net, 0.05)

dim(ginieOutput)

ginieOutput[1:10,]
