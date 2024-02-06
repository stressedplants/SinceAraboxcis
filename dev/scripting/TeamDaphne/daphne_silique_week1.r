#install.packages('Matrix')
library(Matrix)

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#alternative filename='dev/utilities/dataprocessingHelperFunctions.R'
source('dev/utilities/dataprocessingHelperFunctions.R')

a=load('data/GSE226097_silique_230221.RData')  ########CHANGE THIS TO YOUR FILE!!!

araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)

plot(table(sort(clust)), xlab='cluster name', ylab='number of cells', main='Silique')

tfSubs=tfs[which(tfs %in% rownames(gbox))]
length(tfSubs)

thresh=0.01
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))

gbox_filtered=gbox[,includeCells]

dim(gbox_filtered)

numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)

library(umap)

#Let us visualise it using UMAP
gbox.umap <- umap(gbox_filtered)

#Do the cell type clusters group together if we only look at G-box related genes?
colours=rainbow(length(unique(clust))) #YOU MIGHT NEED MORE COLOURS
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='UMAP Silique', xlab='UMAP Component 1', ylab='UMAP Component 2')

pca=prcomp(gbox_filtered, scale.=T, rank.=5)
gbox.pca.umap <- umap(pca$x)

colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]], pch=20, main='PCA UMAP Silique', xlab='UMAP Component 1', ylab='UMAP Component 2')

pdf('figures/Silique_PCAUMAP_by_Cluster.pdf', width=12, height=12)
par(mfrow=c(4,4))
sapply(c(1:length(unique(clust))), function(i){
  colours_specific=rep(rgb(0.1,0.1,0.1,0.1), length(unique(clust)))
  colours_specific[i]='red'
  plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours_specific[clust[includeCells]], pch=20, main=paste('PCA UMAP Silique Cluster', unique(clust)[i]), xlab='UMAP Component 1', ylab='UMAP Component 2')
})
dev.off()


