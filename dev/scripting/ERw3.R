#homework week3 

#install.packages('Matrix')
library(Matrix)
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
source('dev/utilities/dataprocessingHelperFunctions.R')

#reload data 
a=load('data/GSE226097_seedling_3d_230221.RData')
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

#average gene expression of each gene in each cluster 
#make clust not be factors
clustAsNumbers=as.numeric(paste(clust))

geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust))-1), function(j){
    ids=which(clustAsNumbers==j)
    mean(i[ids])
  })
})

colnames(geneExpByCluster)=rownames(gbox)

dim(geneExpByCluster)

#visualizing the outcome of the last analysis as a heatmap
#install.packages('pheatmap')
library('pheatmap')
pheatmap(geneExpByCluster, scale='column')
#cluster labels 
clustLabs=read.table('data/clusterLabels.txt', header=T, sep='\t')
unique(clustLabs[,'Organ'])
organ='Seedlings_3d' 
simpleNames=clustLabs[which(clustLabs[,'Organ']==organ), "Cell.type.suggested"]
print(simpleNames)
rownames(geneExpByCluster)=simpleNames
#redo heatmap with cluster labels 
pheatmap(geneExpByCluster, scale='column')

#save underlying table as a file 
write.table(geneExpByCluster, 'data/Seedling3d_avgExpressionByCluster.txt')

#focussing on transcription factors 
#focussing on transcription factors that bind tp perfect G-boxes. 
#see whether these TFs tend to occur in certain cell types or whether they are more evenly distributed.
tfs=unique(araboxcis[,1])
tfSubs=tfs[which(tfs %in% rownames(gbox))]
pheatmap(geneExpByCluster[,tfSubs], scale='column')

#correlation between expression in TFs and potential downstream targets 
corMat=sapply(tfSubs, function(tf){
  apply(geneExpByCluster, 2, function(gene){
    cor(geneExpByCluster[,tf], gene, method='spearman')
  })
})

dim(corMat)

pheatmap(corMat)

#correlation between TF expression and downstream gene expression in individual cells.
id=which(corMat>0.8 & corMat!=1, arr.ind = TRUE)
dim(id)
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

hist(log(tVal[5,]), main='Seedling3d', xlab='log odds ratio')

#compare between data sets 
#assemble lists of genes that have positive correlations on a cell type level AND a single cell level
thresh=exp(1)
idDoublePositive=which(tVal[5,]>thresh)
doublePositive=id[idDoublePositive,]
doublePositive[,1]=rownames(corMat)[doublePositive[,1]]
doublePositive[,2]=colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive=cbind(doublePositive[,2], doublePositive[,1]) #so TF comes before target 
doublePositive=cbind(doublePositive, tVal[5,idDoublePositive])
#save file
write.table(doublePositive, file='data/Seedling3d_doublePositives.txt', sep='\t', row.names=F) 

#assemble a set of regulatory pairs that were positive in the first analysis and negative in the second analysis.
thresh=1
idSimpson=which(tVal[5,]<thresh)
Simpson=id[idSimpson,]

Simpson[,1]=rownames(corMat)[Simpson[,1]]
Simpson[,2]=colnames(corMat)[as.numeric(Simpson[,2])]
Simpson=cbind(Simpson[,2], Simpson[,1]) #so TF comes before target
Simpson=cbind(Simpson, tVal[idSimpson])
#save file
write.table(Simpson, file='data/Seedling3d_SimpsonPairs.txt', sep='\t', row.names=F)

#making a network from scratch 
#BiocManager::install("GENIE3")
library(GENIE3)
net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=5) #use more than 5 trees (1000)

save(net, file='Seedling3d_network_nTree_5.RData')

ginieOutput=convertToAdjacency(net, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]