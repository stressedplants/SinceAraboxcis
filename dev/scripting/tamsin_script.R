#Week 1 tasks
#--------------------------------------------------

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


#Week 3 tasks
#--------------------------------------------------

#Seurat (Katie's code)
# - Installing Seurat package 
install.packages('Seurat')
library(Seurat)
# - Setting source 
source('dev/utilities/SeuratUMAP.r')
# - Creating Seurat UMAP
SeuratUMAP_function(gbox)

#Finding average gene expression of each gene across each cluster 
# - Converting clust into numeric format (not factors)
clustAsNumbers=as.numeric(paste(clust))
# - Making new matrix with mean expression values for each unique cluster
geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust))-1), function(j){
    ids=which(clustAsNumbers==j)
    mean(i[ids])})})
# - Setting column names 
colnames(geneExpByCluster) = rownames(gbox)
# - Displaying dimensions
dim(geneExpByCluster)

#Visualising data as heatmaps 
# - Installing pheatmap package 
install.packages('pheatmap')
library('pheatmap')
# - Creating heatmap of average gene expression 
pheatmap(geneExpByCluster, scale='column')
# - Reading in the table with cell types 
clustLabs=read.table('data/clusterLabels.txt', header=T, sep='\t')
# - Printing unique names 
unique(clustLabs[,'Organ'])
# - Setting organ name to Seedlings_6d
organ = 'Seedlings_6d'
# - Selecting values from cell.type.suggested column
simpleNames=clustLabs[which(clustLabs[,'Organ']==organ), "Cell.type.suggested"]
# - Printing suggested cell types names 
print(simpleNames)
# - Setting suggested cell types as row names 
rownames(geneExpByCluster)=simpleNames
# - Redoing the heatmap with row names 
pheatmap(geneExpByCluster, scale='column', main='Seedling 6d')
# - Saving the underlying table as a file 
write.table(geneExpByCluster, 'data/Seedlingd6_avgExpressionByCluster.txt')

#Focusing on transcription factors 
# - Storing transcription factor values as variable
tfs = unique(araboxcis[,1])
# - Selecting transcription factors also present in gbox
tfSubs = tfs[which(tfs %in% rownames(gbox))]
# - Creating heatmap with transcription factors 
pheatmap(geneExpByCluster[,tfSubs], scale='column', main='Seedling 6d, TFs')

#Correlation between expression in TFs and potential downstream targets
# - Applying sapply function to calculate Spearman's correlation
# - Storing the results in a matrix 
corMat = sapply(tfSubs, function(tf){
  apply(geneExpByCluster, 2, function(gene){
    cor(geneExpByCluster[,tf], gene, method='spearman')})})
# - Displaying dimensions of the correlation matrix 
dim(corMat)
# - Creating heatmap of correlation matrix
pheatmap(corMat, main='Seedling 6d, Correlation matrix')

#Looking at correlation at the individual cell level 
# Correlation between expression in TFs and 
#  potential downstream targets in individual cells 
# - Creating new matrix with pairs of indices in corMat that have 
#    correlation coefficient > 0.8 and not equal to 1
id = which(corMat>0.8 & corMat!=1, arr.ind = TRUE)
# - Displaying dimensions of correlations > 0.8 
dim(id)
# - Calculating the log-odds ratio using apply function 
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
# - Creating the histogram with appropriate title and axes labels
hist(log(tVal[5,]), main='Seedling 6d', xlab='log odds ratio')

#Comparing between datasets 
#Double positive 
# - Assembling lists of genes that have positive correlations on a cell type level
#     and a single cell level
thresh = exp(1)
idDoublePositive=which(tVal[5,]>thresh)
doublePositive = id[idDoublePositive,]
doublePositive[,1] = rownames(corMat)[doublePositive[,1]]
doublePositive[,2] = colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive=cbind(doublePositive[,2], doublePositive[,1]) #so TF comes before target
doublePositive=cbind(doublePositive, tVal[5,idDoublePositive])
# - Saving double positives as a text file 
write.table(doublePositive, file='data/Seedling6d_doublePositives.txt', sep='\t', row.names=F)
#Simpson pairs 
# - Assembling a set of regulatory pairs that were positive in first analysis
#     and negative in the second analysis (negative regulators?)
thresh=1
idSimpson=which(tVal[5,]<thresh)
Simpson=id[idSimpson,]
Simpson[,1]=rownames(corMat)[Simpson[,1]]
Simpson[,2]=colnames(corMat)[as.numeric(Simpson[,2])]
Simpson=cbind(Simpson[,2], Simpson[,1]) #so TF comes before target
Simpson=cbind(Simpson, tVal[idSimpson])
# - Saving simpson pairs as a text file
write.table(Simpson, file='data/Seedling6d_SimpsonPairs.txt', sep='\t', row.names=F)

#Network analysis from scratch 
# - Importing the Genie3 library 
library(GENIE3)
# - Applying the network algorithm with 5 trees 
net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=5)
save(net, file='seedling6d_network_nTree_5.RData')
load('seedling6d_network_nTree_5.RData')
# - Converting to table format 
ginieOutput=convertToAdjacency(net, 0.05)
# - Displaying dimensions 
dim(ginieOutput)
# - Displaying the first 10 edges in the network
ginieOutput[1:10,]

# - Applying the network algorithm with 50 trees 
net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=50)
save(net, file='seedling6d_network_nTree_50.RData')
load('seedling6d_network_nTree_50.RData')
# - Converting to table format 
ginieOutput=convertToAdjacency(net, 0.05)
# - Displaying dimensions 
dim(ginieOutput)
# - Displaying the first 10 edges in the network
ginieOutput[1:10,]


#Week 5 tasks
#--------------------------------------------------

#Load the network 
# - Loading in the 50 tree network
a = load('seedling6d_network_nTree_50.RData')
# - Accessing the network through GENIE3
newNet = GENIE3::getLinkList(net)

#Finding overlaps between single cell gene network and Araboxcis network
# - Finding unique genes in the new network
genesInNet = unique(c(newNet[,1], newNet[,2]))
# - Filtering the AraBOXcis network to only contain genes that are in the new network
araboxcisFiltered = araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]
# - Extracting the top edges in your new network
newNetTopEdges = newNet[1:length(araboxcisFiltered[,1]),]
# - Reformatting edges to make comparable size to Araboxcis network
edgesNew = paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld = paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep='_')

#Venn Diagram:
# - The overlap
length(which(edgesNew %in% edgesOld))
# - In new network only
length(which(! (edgesNew %in% edgesOld)))
# - In old network only
length(which(! (edgesOld %in% edgesNew)))

#Finding important genes in the network 
# - Exploring which TFs have the most edges/degree
tfsNew = table(newNetTopEdges[,1])
tfsOld = table(araboxcisFiltered[,1])[names(tfsNew)]
# - Plotting histogram of degrees
hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs', main='Seedling 6d, 50 trees')
# - Exploring if the same TFs have high degrees in AraBOXcis and the new network:
plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis', main='Seedling 6d, 50 trees')
# - Displaying 20 TFs with highest degrees  
sort(tfsNew, decreasing=TRUE)[1:20]

#Calculating metrics of TF importance: node betweenness, centrality, hub
# - Loading required libraries
library(igraph)
install.packages('network')
library(network)
library(pheatmap)
# 1) Node betweenness
# - Calculating node betweenness 
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))
node_betweenness_all <- betweenness(simple_network)
node_betweenness=node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing=TRUE)[1:20]
# - Plotting the node betweenness 
plot(sort(node_betweenness), main='Node betweenness')
# 2) Node centrality 
# - Calculating node centrality 
node_centrality_all <- alpha_centrality(simple_network)
node_centrality=node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing=TRUE)[1:20]
# - Plotting node centrality 
plot(sort(node_centrality), main='Node centrality')
# 3) Node hub 
# - Calculating node hub 
node_hub_all <- hub_score(simple_network)$vector
node_hub=node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing=TRUE)[1:20]
# - Plotting node hub
plot(sort(node_hub), main='Node hub')
# - Plotting node betweenness against node centrality 
plot(node_betweenness_all, node_centrality_all)
# - Plotting node hub against node centrality 
plot(node_hub_all, node_centrality_all)
# - Plotting node hub against node betweenness 
plot(node_hub_all, node_betweenness_all)

#Associations between GO terms in the network 
# - Loading data
a = load('C:/Users/tammy/Downloads/functionalData.RData')
source('C:/Users/tammy/OneDrive/Documents/GitHub/SinceAraboxcis/dev/utilities/dataprocessingHelperFunctions.R')
# - Using PAFway to explore GO term associations 
pafwayOut = pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut) = colnames(pafwayOut)
# - Filter to only include rows and columns with at least one significant factor:
atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)
atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)
pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]
# - Creating heatmap of GO terms 
pheatmap(pafwayInterestingOnly)
# - Creating logged heatmap (most significant associations)
pheatmap(log(pafwayInterestingOnly, 10))
