# Loading the data and packages ---------------------------------

install.packages("Matrix")
install.packages("umap")
install.packages("progress")
install.packages("Seurat")
install.packages(c("Rcpp", "ggplot2", "irlba", "Matrix", "splines", "RColorBrewer", 
                   "dplyr", "bipartite", "pracma", "Rtsne", "R.utils"))
BiocManager::install("GENIE3")
library(GENIE3)

library(Seurat)

library(Matrix)
library(progress)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(umap)


source('dev/utilities/dataprocessingHelperFunctions.R') 
#Load the single cell data for seed 12
data = load('data/GSE226097_seedling_12d_230221.RData')
#loading in the og AraBoxcis network 
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)


# What is contained in the data set?  ---------------------------
#Making a frequency plot
#Create the barplot
frequencyplot_clustering_seedling_d12 <- plot(table(sort(clust)),
                                        xlab = 'Cluster Name',
                                        ylab = 'Number of Cells',
                                        main = 'Seedling - Day 12')


# Araboxcis -----------------------------------------------------

#getting a list of tfs
tfs = unique(araboxcis[,1])
tfSubs = tfs[which(tfs %in% rownames(gbox))]
length(tfSubs) #153

#filtering genes and cells
dim(gbox)
thresh = 0.01
numberGenesPerCell = apply(gbox, 2, function(i){length(which(i > 0))})
includeCells = which(numberGenesPerCell > (thresh*dim(gbox)[1]))
gbox_filtered = gbox[,includeCells]
dim(gbox_filtered)

# Vis using UMAP ------------------------------------------------

gbox.umap <- umap(gbox_filtered, n_neighbors = 10)
png("dev/figures/umap_seedling_d12_k10.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], 
     gbox.umap$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'UMAP Seedling d12, k = 10', 
     xlab = 'UMAP Component 1',
     ylab = 'UMAP Component 2',
     xlim = c(-10, 10),
     ylim = c(-7.5, 7.5))
dev.off()

pca = prcomp(gbox_filtered, scale. = T, Rank = 5)

gbox.pca.umap_10 <- umap(pca$x, n_neighbors = 10)
png("dev/figures/umap_PCA_seedling_d12_rud.png", width = 800, height = 600, units = "px")
colours = rainbow(length(unique(clust)))
plot(gbox.pca.umap_10$layout[,1],
     gbox.pca.umap_10$layout[,2], 
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'PCA UMAP Seedling d12, rank undefinded', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')
dev.off()


#Katies code ---------------------------------------

source('dev/utilities/SeuratUMAP.r')
SeuratUMAP_function(gbox_filtered)


# Week 3 --------------------------------------------------------

#alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
source('dev/utilities/dataprocessingHelperFunctions.R')


clustAsNumbers = as.numeric(paste(clust))

geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust)) - 1), function(j){
    ids = which(clustAsNumbers == j)
    mean(i[ids])
  })
})

colnames(geneExpByCluster) = rownames(gbox)

dim(geneExpByCluster)

#install.packages('pheatmap')
library('pheatmap')

pheatmap(geneExpByCluster, scale = 'column')

#By scaling the data in this manner, the pheatmap() function ensures
#that the heatmap representation is not biased towards variables with
#larger magnitudes or variances, allowing for more meaningful
#comparisons between variables or across samples.

clustLabs = read.table('data/clusterLabels.txt',
                       header = T, sep = '\t')

unique(clustLabs[,'Organ'])

organ = 'Seedlings_12d' #### 
simpleNames = clustLabs[which(clustLabs[,'Organ'] == organ),
                        "Cell.type.suggested"]
print(simpleNames)

rownames(geneExpByCluster) = simpleNames

pheatmap(geneExpByCluster, scale = 'column')

#saving the underlying table as a table. 
write.table(geneExpByCluster, 
            'data/Seedling_d12_avgExpressionByCluster.txt')

#Focusing on transcription factor 

tfs = unique(araboxcis[,1])
tfSubs = tfs[which(tfs %in% rownames(gbox))]
pheatmap(geneExpByCluster[,tfSubs], scale = 'column', 
         main = "Heatmap by clustered cell type, Seedlings d12")

#finding correlation between ecery TF and every potential target
corMat = sapply(tfSubs, function(tf){
  apply(geneExpByCluster, 2, function(gene){
    cor(geneExpByCluster[,tf], gene, method = 'spearman')
  })
})

dim(corMat)

# Plot heatmap with axis labels
pheatmap(corMat, 
         main = "Heatmap, TFs and gene correlations, Seedlings d12")

id = which(corMat > 0.8 & corMat != 1, arr.ind = TRUE)
dim(id)

dim(corMat)

#Making a log-odds histogram. 
tVal = apply(id, 1, function(i){
  row = rownames(corMat)[i[1]]
  col = colnames(corMat)[i[2]]
  
  inBoth = length(which(gbox[row,] > 0 & gbox[col,] > 0))
  inNone = length(which(gbox[row,] == 0 & gbox[col,] == 0))
  inFirst = length(which(gbox[row,] > 0 & gbox[col,] == 0))
  inSecond = length(which(gbox[row,] == 0 & gbox[col,] > 0))
  matTemp = matrix(c(inBoth, inFirst, inSecond, inNone), ncol = 2)
  c(inBoth, inFirst, inSecond, inNone, (inBoth*inNone)/(inFirst*inSecond))
})

hist(log(tVal[5,]), main = 'Seedling D12', xlab = 'log odds ratio')

#compare between data sets
thresh = exp(1)
idDoublePositive = which(tVal[5,] > thresh)
doublePositive = id[idDoublePositive,]
doublePositive[,1] = rownames(corMat)[doublePositive[,1]]
doublePositive[,2] = colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive = cbind(doublePositive[,2], doublePositive[,1]) #so TF comes before target
doublePositive = cbind(doublePositive, tVal[5,idDoublePositive])

#save file
write.table(doublePositive, 
            file = 'data/Seedling_d12_doublePositives.txt', 
            sep = '\t', row.names = F)

#data set where regultory paurs that were positive in the first 
# analysis and negartive in the second (negative regulators?)
thresh = 1
idSimpson = which(tVal[5,] < thresh)
Simpson = id[idSimpson,]

Simpson[,1] = rownames(corMat)[Simpson[,1]]
Simpson[,2] = colnames(corMat)[as.numeric(Simpson[,2])]
Simpson = cbind(Simpson[,2], Simpson[,1]) #so TF comes before target
Simpson = cbind(Simpson, tVal[idSimpson])

#save file
write.table(Simpson, file = 'data/Seedling_d12_SimpsonPairs.txt', 
            sep = '\t', row.names = F)

#Network
net_1000 = GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees = 1000) 
save(net, file = 'seedling_d12_network_nTree_1000.RData')

save(net, file = 'seedling_d12_network_nTree_5.RData')

ginieOutput = convertToAdjacency(net, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]
