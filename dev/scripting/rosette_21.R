#install.packages('Matrix')
library(Matrix)

#if (!require("BiocManager", quietly = TRUE))
library("BiocManager")

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
#Priya helped me with this.gbox_filtered <- as.matrix(gbox_filtered)
gbox_matrix <- as.matrix(gbox_filtered)

#Specify parameters

tsne_result <- Rtsne(gbox_filtered, prerplexity = 50, theta = 0.5, dims =2)
tsne_result <- Rtsne(gbox_matrix, dims = 2, perplexity = 45)

#Plot T-sne
plot(tsne_result$Y, col=colours[clust[includeCells]], pch = 20, main = "t-SNE plot Rosette_21d")


#Perform SeuratUMAP
source('dev/utilities/SeuratUMAP.r')
SeuratUMAP_function(gbox_filtered)

#Find the average gene expression of each gene across each cluster
clustAsNumbers=as.numeric(paste(clust))

geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust))-1), function(j){
    ids=which(clustAsNumbers==j)
    mean(i[ids])
  })
})

colnames(geneExpByCluster)=rownames(gbox)

dim(geneExpByCluster)

#Answer: 18 2109

#Load package
library('pheatmap')

#Perform heatmap
pheatmap(geneExpByCluster, scale='column')

#Read in table of what cell type each cluster is
clustLabs=read.table('data/clusterLabels.txt', header=T, sep='\t')

#Print out all uniques names they assign samples:
unique(clustLabs[,'Organ'])

#Indicate organ is Rosette
organ='Rosette 21d'
simpleNames=clustLabs[which(clustLabs[,'Organ']==organ), "Cell.type.suggested"]
print(simpleNames)

#assign rownames
rownames(geneExpByCluster)=simpleNames

#Perform heatmap again
pheatmap(geneExpByCluster, scale='column', main="Rosette 21D")

#Save file
write.table(geneExpByCluster, 'data/Rosette21d_avgExpressionByCluster.txt')

#Perform heatmap but this time the analysis looks at trasnciption factors
tfs=unique(araboxcis[,1])
tfSubs=tfs[which(tfs %in% rownames(gbox))]
pheatmap(geneExpByCluster[,tfSubs], scale='column', main= "Rosette 21D")

#Find correlation between every TF and every potential downstream target across different cell types.
corMat=sapply(tfSubs, function(tf){
  apply(geneExpByCluster, 2, function(gene){
    cor(geneExpByCluster[,tf], gene, method='spearman')
  })
})

dim(corMat)
#Answer: 2109  145

#Create heatmap
pheatmap(corMat)

#Calcualte correlation the TFs that are highly correlated with downstream genes
id=which(corMat>0.8 & corMat!=1, arr.ind = TRUE)
dim(id)
#Answer= 1331 2

#Check corMat and Id values
dim(corMat)
#Answer: 2109 2
dim(id)
#Answer: 1331 2
 
#Perform Histogram
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

hist(log(tVal[5,]), main='Rosette 21D', xlab='log odds ratio')



#assemble lists of genes that have positive correlations on a cell type level AND a single cell level.
thresh=exp(1)
idDoublePositive=which(tVal[5,]>thresh)
doublePositive=id[idDoublePositive,]
doublePositive[,1]=rownames(corMat)[doublePositive[,1]]
doublePositive[,2]=colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive=cbind(doublePositive[,2], doublePositive[,1]) #so TF comes before target
doublePositive=cbind(doublePositive, tVal[5,idDoublePositive])

#save file
write.table(doublePositive, file='data/Rosette21D_doublePositives.txt', sep='\t', row.names=F) ##CHANGE FILE NAME

#Assemble a set of regulatory pairs that were positive in the first analysis and negative in the second analysis.
thresh=1
idSimpson=which(tVal[5,]<thresh)
Simpson=id[idSimpson,]

Simpson[,1]=rownames(corMat)[Simpson[,1]]
Simpson[,2]=colnames(corMat)[as.numeric(Simpson[,2])]
Simpson=cbind(Simpson[,2], Simpson[,1]) #so TF comes before target
Simpson=cbind(Simpson, tVal[idSimpson])

#save file
write.table(Simpson, file='data/Rosette21D_SimpsonPairs.txt', sep='\t', row.names=F)

#Load package
library(GENIE3)

# Start timing
start_time <- Sys.time()

# Code to perform genie3 
net <- GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees = 50) # Adjust nTrees as needed

# Calculate elapsed time
elapsed_time <- Sys.time() - start_time

# Print elapsed time
print(elapsed_time)
#Anser: 10.32947 hours

# Save the object
save(net, file = 'Rosette21D_network_nTree_50.RData')

#Convert the 50 trees into a table 
ginieOutput=convertToAdjacency(net, 0.05)

#Check dimension
dim(ginieOutput)
#Answer: 257   3

#Table output
ginieOutput[1:10,]
