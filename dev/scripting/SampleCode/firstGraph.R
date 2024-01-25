install.packages('pheatmap')
install.packages('umap')
library(pheatmap)
library(umap)

#Load the single cell data for flowers.  
load('data/GSE226097_flower_230221.RData')

#gbox is a Matrix containing the expression for all 2077 genes from the AraBOXcis network for each cell
dim(gbox)
rownames(gbox)

#clust is a vector containing the 'cluster' of the cell.  
#Cells that are in the same cluster are probably the same type (like two epidermis cells...)
length(clust)
#How many cells are in each cluster?
table(clust)
plot(table(clust))

#Let us make a smaller Matrix that just contains cells from one cluster
cluster=12
subMat=gbox[,which(clust==cluster)]
dim(subMat)

#remove genes with 0 expression
subMat=subMat[which(rowSums(subMat)!=0),]

#visualise heatmap
pheatmap(subMat, scale='row')

#Let us visualise it using UMAP
gbox.umap <- umap(gbox)
plot(gbox.umap$layout[,1], gbox.umap$layout[,2])

#Do the cell type clusters group together if we only look at G-box related genes?
colours=rainbow(12)
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust], pch=20)

#What if we change the parameters of UMAP?
gbox.umap <- umap(gbox, n_neighbors=10)
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust], pch=20)

