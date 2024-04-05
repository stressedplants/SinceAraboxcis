#Load in the three networks, renaming each one after it is loaded,
#as they all were called net. 
a = load("data/seedling_d12_network_nTree_50.RData")
net_d12_50 = net_50
b = load("data/seedling6d_network_nTree_50.RData")
net_d06_50 = net
c = load("data/Seedling3d_network_nTree_50.RData")
net_d03_50 = net


#converts the weight matrix 
d12Net = GENIE3::getLinkList(net_d12_50)
d06Net = GENIE3::getLinkList(net_d06_50)
d03Net = GENIE3::getLinkList(net_d03_50)

#(the gene network) into a sorted list of regulatory links. 
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)

#gets the set of unique genes in each network
uniquegenes_d12Net = unique(c(d12Net[,1], d12Net[,2]))
uniquegenes_d06Net = unique(c(d06Net[,1], d06Net[,2]))
uniquegenes_d03Net = unique(c(d03Net[,1], d03Net[,2]))

#Findiing all the common genes
common_genes <- Reduce(intersect, 
                       list(uniquegenes_d12Net, uniquegenes_d06Net, 
                            uniquegenes_d03Net))

#Creating filtered data sets where the rows only have 
#genes that are common to all data sets
d12Net_filtered = d12Net[which(d12Net$regulatoryGene %in% common_genes
                              & d12Net$targetGene %in% common_genes),]

d06Net_filtered = d12Net[which(d06Net$regulatoryGene %in% common_genes
                               & d06Net$targetGene %in% common_genes),]

d03Net_filtered = d12Net[which(d03Net$regulatoryGene %in% common_genes
                               & d03Net$targetGene %in% common_genes),]

#reformatting so edges can be compared(this turn the data into 1 column)
edges_d12 = paste(d12Net_filtered[,1], d12Net_filtered[,2], sep = '_')
edges_d06 = paste(d06Net_filtered[,1], d06Net_filtered[,2], sep = '_')
edges_d03 = paste(d03Net_filtered[,1], d03Net_filtered[,2], sep = '_')

#finding the edges in all dX using intsect
length(intersect(intersect(edges_d12, edges_d06), edges_d03))
common_edges <- intersect(intersect(edges_d12, edges_d06), edges_d03)

#finding the edges in all two of the three data sets using intsect
d12_d06_common_edges <- intersect(edges_d12, edges_d06)
#removing the edges that are in the common_edges data set
d12_d06_common_edges <- setdiff(d12_d06_common_edges, common_edges)
length(d12_d06_common_edges)

#finding the edges in all two of the three data sets using intsect
d12_d03_common_edges <- intersect(edges_d12, edges_d03)
#removing the edges that are in the common_edges data set
d12_d03_common_edges <- setdiff(d12_d03_common_edges, common_edges)
length(d12_d03_common_edges)

#finding the edges in all two of the three data sets using intsect
d06_d03_common_edges <- intersect(edges_d06, edges_d03)
#removing the edges that are in the common_edges data set
d06_d03_common_edges <- setdiff(d06_d03_common_edges, common_edges)
length(d06_d03_common_edges)

#finding the edges only in each dataset
length(setdiff(setdiff(setdiff(edges_d12, common_edges), d12_d03_common_edges),d12_d06_common_edges )) 

length(setdiff(setdiff(setdiff(edges_d06, common_edges), d06_d03_common_edges), d12_d06_common_edges))
length(setdiff(setdiff(setdiff(edges_d03, common_edges), d06_d03_common_edges), d12_d03_common_edges))
