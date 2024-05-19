#This section of code is used to generate a plot showing the number of cells in each cluster in the snRNA-seq data of 30-day-old rosettes of Arabidopsis

library(Matrix)

#Loads file which contains functions to help analyse the code
source('dev/utilities/dataprocessingHelperFunctions.R')


#Loads the snRNA-seq data of 30-day-old rosettes of Arabidopsis
a = load('data/GSE226097_rosette_30d_230221.RData')


#Loads the original AraBOXcis network based on the bulk RNA-seq data in 7-day-old seedlings
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)


#This will show the variables stored in the single cell data. This data has two variables gbox and clust
print(a)


#Tells the number of rows (these are the G-box related genes) and the number of columns (these are the individual cells) in gbox. The output is 2098, 37817 so there are 2098 genes and 37817 cells.
dim(gbox)


#Shows the length of the variable cluster, 37817 which is the number of cells. Each cell in gbox has a cluster designation.
length(clust)


#Converts gbox into a Matrix i.e. a normal table in R
as.matrix(gbox[1:50, 1:3])


#Function table counts the number of cells in each cluster
table(clust)

#This plots the number of cells in each cluster
plot(table(sort(clust)), xlab='Cluster Designation', ylab='No. of Cells', main='No. of cells in each cluster in gbox')
#---------------------------------------------------------------------
#This section is used to calculate how many unique G-box binding TFs there are in gbox


#Shows the dimensions of the bulk RNA-seq data used for the AraBOXcis network
dim(araboxcis)
#There are 5000 rows and 3 columns

#Shows the first 4 rows of araboxcis
araboxcis[1:4,]

#column 1 is the TF
#column 2 is the target gene
#column 3 is the score of the regulatory edge


#Makes a histogram of the scores of the regulatory edges in araboxcis
hist(araboxcis[,3])


#Stores a list of the unique TFs in araboxcis in tfs
tfs=unique(araboxcis[,1])


#tfSubs is filtered to only include the TFs that are present in both gbox and araboxcis
tfSubs=tfs[which(tfs %in% rownames(gbox))]


#Finds out how many TFs (144) there are in gbox.
length(tfSubs)

#---------------------------------------------------------------------
#This section of code was used to filter the genes and cells in gbox with low values


#Sets the threshold to 1%
thresh=0.01


#This function counts the no. of genes in each cell where the gene expression  bigger than 0
numberGenesPerCell=apply(gbox, 2, function(i){length(which(i>0))})


#This filters gbox to remove the cells that have less than 1% of the genes expressed at all
includeCells=which(numberGenesPerCell>(thresh*dim(gbox)[1]))


#The filtered gbox is stored in gbox_filtered
gbox_filtered=gbox[,includeCells]


#Looks at the dimensions of g+box filtered. There are 2098 genes (rows) and 37816 (cells). This means that only 1 cell was removed.
dim(gbox_filtered)


#This filters gbox_filtered to remove the genes expressed in less than 1% of the cells
numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))


gbox_filtered=gbox_filtered[includeGenes,]


#Shows dimensions of gbox_filtered. After filtering, there are 1520 genes (Rows) and 37816 cells (columns).
dim(gbox_filtered)
#---------------------------------------------------------------------
#This section of code was used to generate the PCA before UMAP plot to visualise the snRNA-seq data
#UMAP

#Load the required packages
install.packages('umap')
library(umap)

#This code will retrive data from clusterLabels.txt which includes the cell type annotation for each cluster
clustLabs=read.table('data/clusterLabels.txt', header=T, sep='\t')
unique(clustLabs[,'Organ'])
organ='Rosette 30d'
simpleNames=clustLabs[which(clustLabs[,'Organ']==organ), "Cell.type.suggested"]


#This shows all the cell type annotations for gbox
print(simpleNames)


#This will plot a PCA then UMAP plot
pca <- prcomp(gbox_filtered, scale. = TRUE, rank. = 10)
gbox.pca.umap <- umap(pca$x)

colours <- c("hotpink", "blue", "green", "green4", "purple", "orange", "yellow", "cyan",
             "magenta", "brown", "pink", "darkblue", "darkorchid4", "red",
             "springgreen", "darkorange", "darkcyan")

colour_mapping <- setNames(colours, simpleNames)

plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], 
     col = colour_mapping[clust[includeCells]], pch = 20, 
     main = 'PCA UMAP Rosette 30 Days', 
     xlab = 'UMAP Component 1', ylab = 'UMAP Component 2')

# Key for the cluster designations in the plot
# 0 hotpink = unannotated
# 1 blue = unnannotated
# 2 green = adaxial epidermal
# 3 green4 = adaxial epidermal
# 4 purple = vascular
# 5 orange = companion cells
# 6 yellow = unannotated
# 7 cyan = Procambium_PP
# 8 magenta = Unannotated
# 9 brown = Phloem
# 10 pink = Unannotated
# 11 darkblue  Trichome
# 12 darkorchid4 = Mesophyll
# 13 red = Xylem
# 14 springgreen = Epidermal
# 15 darkorange = Dividing
# 16 darkcyan = Guard


#This code was used to generate PCA then UMAP plots to look at individual cell types
colours <- c("grey83", "grey83", "grey83", "grey83", "grey83", "grey83", "grey83", "grey83",
             "grey83", "grey83", "grey83", "grey83", "grey83", "grey83", "grey83", "grey83", "darkcyan")

colour_mapping <- setNames(colours, simpleNames)

plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], 
     col = colour_mapping[clust[includeCells]], pch = 20, 
     xlab = 'UMAP Component 1', ylab = 'UMAP Component 2')
#---------------------------------------------------------------------
#This section describes how the GRN was made for the snRNA-seq data of 30-day-old rosettes


#Load required packages. GENIE3 is the network inference package we used.
library(GENIE3)


#The parameters for the network were: The input is the gene expression dataset, gbox, the regulator is the list of G-box binding TFs, and the number of trees in the random forest was 50.
net=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees=50)


save(net, file='../../rosette30d_network_nTree_50.RData')


#Converts R data file to a table. We used the convertToAdjacency function from dataprocessingHelperFunctions.R
ginieOutput=convertToAdjacency(net, 0.05)


dim(ginieOutput)


ginieOutput[1:10,]
#---------------------------------------------------------------
#This section will find the average gene expression of each gene across the different cell types in gbox in a heatmap

#Load required packages
install.packages('pheatmap')
library('pheatmap')


#Converts clust (the cluster designation of each cell) into a numeric vector
clustAsNumbers=as.numeric(paste(clust))


#Calculates the average expression of each gene within each cluster
geneExpByCluster = apply(gbox, 1, function(i){
  sapply(0:(length(unique(clust))-1), function(j){
    ids=which(clustAsNumbers==j)
    mean(i[ids])
  })
})

#Sets column names of the geneExpByCluster to be the same as the row names of gbox 
colnames(geneExpByCluster)=rownames(gbox)

#Shows the no. of rows (clusters = 17) + columns (genes =2098) of the geneExpByCluster matrix
dim(geneExpByCluster)

#Visualising large tables as heatmaps
install.packages('pheatmap')
library('pheatmap')


#Retrieve the file which contains the cell type annotation for each cluster
clustLabs=read.table('data/clusterLabels.txt', header=T, sep='\t')
unique(clustLabs[,'Organ'])
organ='Rosette 30d'
simpleNames=clustLabs[which(clustLabs[,'Organ']==organ), "Cell.type.suggested"]
print(simpleNames)

#Assigns each cluster to a cell type
rownames(geneExpByCluster)=simpleNames


#Gets all the unique TFs in araboxcis
tfs=unique(araboxcis[,1])


#Filters to only include the TFs that are present in both araboxcis and gbox
tfSubs=tfs[which(tfs %in% rownames(gbox))]


#Generates heatmap showing average gene expression of G-box binding TFs across the cell types
pheatmap(geneExpByCluster[,tfSubs], scale='column', main='Rosette 30d')

#---------------------------------------------------------------------
#This section looks at the correlation between the expression of G-box binding TFs and the potential downstream targets


#Generates heatmap showing correlation between the expression of G-box binding TFs and the potential downstream targets using Spearman's correlation
corMat=sapply(tfSubs, function(tf){
  apply(geneExpByCluster, 2, function(gene){
    cor(geneExpByCluster[,tf], gene, method='spearman')
  })
})


dim(corMat) #2098 144


pheatmap(corMat, main='Rosette 30d')

#---------------------------------------------------------------------
#This section of code looks at the correlation between TF expression + downstream target genes in individual cells


#Filters the correlation matrix for where value is >0.8
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

#Generates histogram displaying the log odds ratio of observing the TF and the downstream gene within the same cell
hist(log(tVal[5,]), main='Rosette 30d', xlab='log odds ratio')

#---------------------------------------------------------------------
#This section of code generates list of TF-target gene interactions that show correlation at both a cell type level and an individual level
thresh=exp(1)
idDoublePositive=which(tVal[5,]>thresh)
doublePositive=id[idDoublePositive,]
doublePositive[,1]=rownames(corMat)[doublePositive[,1]]
doublePositive[,2]=colnames(corMat)[as.numeric(doublePositive[,2])]
doublePositive=cbind(doublePositive[,2], doublePositive[,1])
doublePositive=cbind(doublePositive, tVal[5,idDoublePositive])


#save file
write.table(doublePositive, file='data/Rosette30d_doublePositives.txt', sep='\t', row.names=F)


#This section of code generates list of TF-target gene interactions that show opposite correlation at a cell type level and an individual level; showing Simpson's Paradox
thresh=1
idSimpson=which(tVal[5,]<thresh)
Simpson=id[idSimpson,]


Simpson[,1]=rownames(corMat)[Simpson[,1]]
Simpson[,2]=colnames(corMat)[as.numeric(Simpson[,2])]
Simpson=cbind(Simpson[,2], Simpson[,1])
Simpson=cbind(Simpson, tVal[idSimpson])


#save file
write.table(Simpson, file='data/Rosette30d_SimpsonPairs.txt', sep='\t', row.names=F)

#---------------------------------------------------------------------
#This section of code was used to find the overlapping TF-target gene regulatory edges/interactions between the previous araboxcis network and the new gbox network. Also to see if the new network is scale-free.

#Loads the networks
a=load('rosette30d_network_nTree_50.RData')
newNet=GENIE3::getLinkList(net)
araboxcis=read.csv('data/gboxNetwork22C.csv', header=T)


genesInNet=unique(c(newNet[,1], newNet[,2]))


araboxcisFiltered=araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]


newNetTopEdges=newNet[1:length(araboxcisFiltered[,1]),]


edgesNew=paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld=paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep='_')


#Calculates the overlap: 9307
length(which(edgesNew %in% edgesOld)) #overlap 9307


#Calculates the edges unique to gbox = 37425
length(which(! (edgesNew %in% edgesOld))) #new network only 37425


#Calculates the edges unique to araboxcis = 37425
length(which(! (edgesOld %in% edgesNew)))


tfsNew=table(newNetTopEdges[,1])
tfsOld=table(araboxcisFiltered[,1])[names(tfsNew)]


#sPlots histogram to see if the new network is scale-free
hist(as.numeric(tfsNew), main='SinceAraBOXcis', xlab='degree of TFs')

plot(as.numeric(tfsNew), as.numeric(tfsOld), xlab='degree in SinceAraBOXcis', ylab='degree in AraBOXcis')


#Prints out the top 20 edges in gbox
sort(tfsNew, decreasing=TRUE)[1:20]

#--------------------------------------------------------------------
#This section was used to look at the metrics(hub_score and betweenness) of TF importance

#Load required packages
library(igraph)
library(network)


simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))


node_betweenness_all <- betweenness(simple_network)
node_betweenness=node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing=TRUE)[1:20]


plot(sort(node_betweenness))


node_hub_all <- hub_score(simple_network)$vector
node_hub=node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing=TRUE)[1:20]


plot(sort(node_hub))


plot(node_hub_all, node_betweenness_all)
#---------------------------------------------------------------------
#This section found associations between GO terms in the network

library(pheatmap)
a=load('data/functionalData.RData')
source('dev/utilities/dataprocessingHelperFunctions.R')


#The statistical tool pafway calculates the p-value of the association between GO terms
pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut)=colnames(pafwayOut)


atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)


atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)


pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]

#Generates heatmap of the GO terms in columns that act upstream of the Go terms in the rows. The values correspond to the log10 of the p-values.
pheatmap(log(pafwayInterestingOnly, 10))

#---------------------------------------------------------------------
#This section calculates the top 1000 overlapping edges between araboxcis and gbox to be entered into cytoscape


#The value should have 9307 elements, as there are 9307 overlapping edges between gbox and araboxcis
overlap <- edgesNew[edgesNew %in% edgesOld]


overlap_df <- data.frame(overlap = overlap,
                         regulatoryGene = newNetTopEdges$regulator[match(overlap, edgesNew)],
                         targetGene = newNetTopEdges$target[match(overlap, edgesNew)],
                         weight = newNetTopEdges$weight[match(overlap, edgesNew)])


write.table(overlap_df[1:1000, -1], file='overlaptop100030drosetteandaraboxcis.csv', sep=',', row.names=F)

#---------------------------------------------------------------------
#This section calculates the top 1000 overlapping edges between the 21d-old rosette and 30d-old rosette (gbox) to be entered into Cytoscape


# Loads the 21d and 30d rosette datasets
load('rosette30d_network_nTree_50.RData')
net_30d <- GENIE3::getLinkList(net)


load('Rosette21d_network_nTree_50.RData')
net_21d <- GENIE3::getLinkList(net)


#Finds the overlapping edges
overlap21dand30d <- intersect(paste(net_30d[, 1], net_30d[, 2], sep = '_'),
                              paste(net_21d[, 1], net_21d[, 2], sep = '_'))


overlap_21dand30d <- data.frame(overlap = overlap21dand30d,
                         regulator_30d = net_30d$regulator[match(overlap21dand30d, paste(net_30d[, 1], net_30d[, 2], sep = '_'))],
                         target_30d = net_30d$target[match(overlap21dand30d, paste(net_30d[, 1], net_30d[, 2], sep = '_'))],
                         weight_30d = net_30d$weight[match(overlap21dand30d, paste(net_30d[, 1], net_30d[, 2], sep = '_'))],
                         regulator_21d = net_21d$regulator[match(overlap21dand30d, paste(net_21d[, 1], net_21d[, 2], sep = '_'))],
                         target_21d = net_21d$target[match(overlap21dand30d, paste(net_21d[, 1], net_21d[, 2], sep = '_'))],
                         weight_21d = net_21d$weight[match(overlap21dand30d, paste(net_21d[, 1], net_21d[, 2], sep = '_'))])



write.csv(overlap_21dand30d[1:1000, c(2, 3, 4)], file='overlaptop1000of21dand30drosette.csv', row.names = FALSE)

#---------------------------------------------------------------------
#This section calculates the number of edges that overlap between the 21d and 30d rosettes


genesInNet=unique(c(newNet[,1], newNet[,2]))


#Filters the 21d rosette data to include the genes that are also in the 30d network
net_21dFiltered=net_21d[which(net_21d[,1] %in% genesInNet & net_21d[,2] %in% genesInNet),]


#Makes the 21d and 30d rosette networks the same size
newNetTopEdges=newNet[1:length(net_21dFiltered[,1]),]


edgesNew=paste(newNetTopEdges[,1], newNetTopEdges[,2], sep='_')
edgesOld=paste(net_21dFiltered[,1], net_21dFiltered[,2], sep='_')


#Calculates the overlapping edges = 296449
length(which(edgesNew %in% edgesOld))


#Calculates the number of edges In 30d rosette network only = #2707
length(which(! (edgesNew %in% edgesOld)))


#Calculates the number of edges In 21d rosette network only = #2707
length(which(! (edgesOld %in% edgesNew)))