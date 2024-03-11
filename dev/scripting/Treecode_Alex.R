# Loading the data and packages ---------------------------------

#install.packages("Matrix")
#install.packages("umap")
#install.packages("progress")
#install.packages("Seurat")
#install.packages(c("Rcpp", "ggplot2", "irlba", "Matrix", "splines", "RColorBrewer", 
#                   "dplyr", "bipartite", "pracma", "Rtsne", "R.utils"))
#BiocManager::install("GENIE3")
library(GENIE3)
library(Matrix)


source('dev/utilities/dataprocessingHelperFunctions.R') 
#Load the single cell data for seed 12
data = load('data/GSE226097_seedling_12d_230221.RData')
#loading in the og AraBoxcis network 
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)


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


#Network

net_50 = GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees = 50) 
save(net_50, file = 'data/seedling_d12_network_nTree_50.RData')


net_100 = GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees = 100) 
save(net_100, file = 'data/seedling_d12_network_nTree_100.RData')


ginieOutput = convertToAdjacency(net_5, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]
