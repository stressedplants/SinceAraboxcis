# aims of this week: 
## how does the average gene expression differ bt cell types? 
## what is the corr. between the expression of g-box binding tfs and genes found near gbox motifs
## how do the correlations calc. above vary if the data is grouped by cell type vs. ind. cell
library(monocle3)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(umap)
library(dplyr)
library(tidyverse)
library(BiocManager)
library(Matrix)
library(patchwork)
library(R.utils)
library(shiny)
library(SingleCellExperiment)
library(mclust)

library(Matrix)
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")  
source('dev/utilities/dataprocessingHelperFunctions.R')
hw2data = load('data/GSE226097_seed_0d_230221.RData')
araboxcis = read.csv("data/gboxNetwork22C.csv", header = T)

# find out avg gene expression of each gene across each cluster

#make clust not be factors
clustAsNumbers <- as.numeric(paste(clust))

geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust))-1), function(j){
    ids=which(clustAsNumbers==j)
    mean(i[ids])
  })
})

colnames(geneExpByCluster) = rownames(gbox)

dim(geneExpByCluster)

#install the heatmap package, load it, and generate the first heatmap
install.packages('pheatmap')
library('pheatmap')
pheatmap(geneExpByCluster, scale = 'column')

#read in the cluster labels provided by the researchers
clustLabs = read.table('data/clusterLabels.txt', header = T, sep = '\t')

#print out all unique names assigned to samples: 
unique(clustLabs[, 'Organ'])

#add the cluster names for my specific plant organ to the data
organ <- 	
  "Seed_0d"
simpleNames <- clustLabs[which(clustLabs[,'Organ'] == organ), "Cell.type.suggested"]
print(simpleNames)
rownames(geneExpByCluster) <- simpleNames

#redo the heatmap with cluster names

pheatmap(geneExpByCluster, scale = 'column')

#save the data as a file: 

write.table(geneExpByCluster, 'data/Seed_0d_avgExpressionByCluster.txt')

#Focus on TF's: 

tfs <- unique(araboxcis[,1])
tfSubs <- tfs[which(tfs %in% rownames(gbox))]
pheatmap(geneExpByCluster[,tfSubs], scale = 'column')

#find correlation between TF's and downstream targets across cell types 

corMat <- sapply(tfSubs, function(tf){
  apply(geneExpByCluster, 2, function(gene){
    cor(geneExpByCluster[,tf], gene, method='spearman')
  })
})
dim(corMat)
pheatmap(corMat)


#correlation between tf's and downstream targets in individual cells
id=which(corMat>0.8 & corMat!=1, arr.ind = TRUE)
dim(id)
dim(corMat)
tVal <- apply(id, 1, function(i){
  row = rownames(corMat)[i[1]]
  col = colnames(corMat)[i[2]]
  
  inBoth = length(which(gbox[row,]>0 & gbox[col,]>0))
  inNone = length(which(gbox[row,] == 0 & gbox[col,] == 0))
  inFirst = length(which(gbox[row,]>0 & gbox[col,] == 0))
  inSecond = length(which(gbox[row, ] == 0 & gbox[col,]>0))
  matTemp = matrix(c(inBoth, inFirst, inSecond, inNone), ncol = 2)
  c(inBoth, inFirst, inSecond, inNone, (inBoth*inNone)/(inFirst*inSecond))
})

hist(log(tVal[5,]), main = "Seed_0d", xlab = 'log odds ratio')

thresh = exp(1)
idDoublePositive <- which(tVal[5,]>thresh)
doublePositive <- id[idDoublePositive,]
doublePositive[,1] <- rownames(corMat)[doublePositive[,1]]
doublePositive[,2] <- colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive <- cbind(doublePositive[,2], doublePositive[,1])
doublePositive <- cbind(doublePositive, tVal[5, idDoublePositive])

#savefile
write.table(doublePositive, file = 'data/Seed0d_doublePositives.txt', sep = '\t', row.names = F)

thresh = 1
idSimpson <- which(tVal[5,]< thresh)
Simpson <- id[idSimpson,]

Simpson[,1] <- rownames(corMat)[Simpson[,1]]
Simpson[,2] <- colnames(corMat)[as.numeric(Simpson[,2])]
Simpson <- cbind(Simpson[,2], Simpson[,1])
Simpson <- cbind(Simpson, tVal[idSimpson])

#save file
write.table(Simpson, file = 'data/Seed0d_SimpsonPairs.txt', sep = '\t', row.names = F)


#making a network
library(GENIE3)
#BiocManager::install("GENIE3")

net <- GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees = 5)
save(net, file = 'seed0d_network_nTree_5.RData')

net1000 <- GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees = 1000)
save(net1000, file = 'seed0d_network_nTree_1000.RData')

#convert this table 
genieOutput <- convertToAdjacency(net1000, 0.05)
dim(genieOutput)
genieOutput[1:10,]
