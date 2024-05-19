
# Seedling overlap, using Net size of araboxcis overlap ---------

#Load in the three networks, renaming each one after it is loaded,
#as they all were called net. 
a = load("data/seedling_d12_network_nTree_50.RData")
net_d12_50 = net_50
b = load("data/seedling6d_network_nTree_50.RData")
net_d06_50 = net
c = load("data/Seedling3d_network_nTree_50.RData")
net_d03_50 = net
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)

#converts the weight matrix 
d12Net = GENIE3::getLinkList(net_d12_50)
d06Net = GENIE3::getLinkList(net_d06_50)
d03Net = GENIE3::getLinkList(net_d03_50)

#I am finding the unique gene in each network and then the unquie genes in all of them. 
genesInNet_d12 = unique(c(d12Net[,1], d12Net[,2]))
genesInNet_d06 = unique(c(d06Net[,1], d06Net[,2]))
genesInNet_d03 = unique(c(d03Net[,1], d03Net[,2]))

combined_genes <- union(union(genesInNet_d12, genesInNet_d06), genesInNet_d03)

#filter araboxcis for these genes.
araboxcisFiltered = araboxcis[which(araboxcis[,1] %in% combined_genes & araboxcis[,2] %in% combined_genes),]

#Now resize the networks to match araboxcis. 5
d12_NetTopEdges = d12Net[1:13895,]
d06_NetTopEdges = d06Net[1:13895,]
d03_NetTopEdges = d03Net[1:13895,]



#I am making a vectors of [regulator gene name]_[target gene name]
#This will make it easier to find the commmon, and exclusive edges. 
edges_d12 = paste(d12_NetTopEdges[,1], d12_NetTopEdges[,2], sep = '_')
edges_d06 = paste(d06_NetTopEdges[,1], d06_NetTopEdges[,2], sep = '_')
edges_d03 = paste(d03_NetTopEdges[,1], d03_NetTopEdges[,2], sep = '_')
edges_araboxcis = paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep = '_')

# intersect() function is used to find the common elements between two
# or more vectors or sets. It returns a vector containing only the
# elements that are present in all the input vectors. 

#Here I use it to find all the common edges in all three networks
common_edges <- intersect(intersect(intersect(edges_d12, edges_d06), edges_d03), edges_araboxcis) 
common_edges_for_net <- common_edges

#filtering the edges for those only found in araboxcis nextwork 
ara_edges_d12 <- intersect(edges_d12, edges_araboxcis)
ara_edges_d06 <- intersect(edges_d06, edges_araboxcis)
ara_edges_d03 <- intersect(edges_d03, edges_araboxcis)

#print total number of common edges
length(common_edges)
# setdiff() function is used to find the set difference between
#two sets. It returns a vector containing elements that are present in
#the first set but not in the second set.

#Here I used intersect to find the common edges between two networks, and 
#then setdiff to remove the one already in common_edges
exclusive_edges_d12_d06 <- setdiff(intersect(ara_edges_d12, ara_edges_d06), common_edges)
#print total number of common edges
length(exclusive_edges_d12_d06)

#I repeat this for the other intersections 
exclusive_edges_d12_d03 <- setdiff(intersect(ara_edges_d12, ara_edges_d03), common_edges)
length(exclusive_edges_d12_d03)

exclusive_edges_d03_d06 <- setdiff(intersect(ara_edges_d03, ara_edges_d06), common_edges)
length(exclusive_edges_d03_d06)

#I use setdiff twice to find all the unqiue egdes in one networks and then print its value
exclusive_edges_d12 <- setdiff(setdiff(ara_edges_d12, ara_edges_d06), edges_d03)
length(exclusive_edges_d12)
exclusive_edges_d06 <- setdiff(setdiff(ara_edges_d06, ara_edges_d12), edges_d03)
length(exclusive_edges_d06)
exclusive_edges_d03 <- setdiff(setdiff(ara_edges_d03, ara_edges_d06), edges_d12)
length(exclusive_edges_d03)

#I am compiling the data into a table. Each of the following will be a row,

d03_data <- c(length(exclusive_edges_d03), length(exclusive_edges_d03_d06),
              length(exclusive_edges_d12_d03), length(common_edges))
d06_data <- c(length(exclusive_edges_d03_d06), length(exclusive_edges_d06),
              length(exclusive_edges_d12_d06), length(common_edges))
d12_data <- c(length(exclusive_edges_d12_d03), length(exclusive_edges_d12_d06),
              length(exclusive_edges_d12), length(common_edges))

Venn_matrix <- matrix(c(d03_data, d06_data, d12_data), nrow = 3,
                      ncol = 4, byrow = TRUE, 
                           dimnames = list(c("d03", "d06", "d12"), 
                                           c("shared_d03", "shared_d06",
                                             "shared_d12", "common")))
# NOTE: where d03 meets shared_d03 it is the number of unique edges in that network. 

#add new column of sum of rows, I can see that this matches the 
# original length of each network
Venn_matrix <- cbind(Venn_matrix, rowSums(Venn_matrix))

#rename new column
colnames(Venn_matrix) <- c(colnames(Venn_matrix)[1:4], "total")

# Print the updated matrix
print(Venn_matrix)
write.csv(Venn_matrix, file = "venndiagram_ara.csv", row.names = FALSE)

#making a list of the edges and their weight that are common do d03,d06, d12 and araboxcis

commonNet_genes <- strsplit(common_edges_for_net , "_")
regulatoryGene <- sapply(commonNet_genes, "[[", 1)
targetGene <- sapply(commonNet_genes, "[[", 2)
commonNet <- data.frame(regulatoryGene, targetGene)

write.csv(commonNet, "data/commonNet.csv", row.names = FALSE)
