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

#I am making a vectors of [regulator gene name]_[target gene name]
#This will make it easier to find the commmon, and exclusive edges. 
edges_d12 = paste(d12Net[,1], d12Net[,2], sep = '_')
edges_d06 = paste(d06Net[,1], d06Net[,2], sep = '_')
edges_d03 = paste(d03Net[,1], d03Net[,2], sep = '_')

# intersect() function is used to find the common elements between two
# or more vectors or sets. It returns a vector containing only the
# elements that are present in all the input vectors. 

#Here I use it to find all the common edges in all three networks
common_edges <- intersect(intersect(edges_d12, edges_d06), edges_d03)
#print total number of common edges
length(common_edges)

# setdiff() function is used to find the set difference between
#two sets. It returns a vector containing elements that are present in
#the first set but not in the second set.

#Here I used intersect to find the common edges between two networks, and 
#then setdiff to remove the one already in common_edges
exclusive_edges_d12_d06 <- setdiff(intersect(edges_d12, edges_d06), common_edges)
#print total number of common edges
length(exclusive_edges_d12_d06)

#I repeat this for the other intersections 
exclusive_edges_d12_d03 <- setdiff(intersect(edges_d12, edges_d03), common_edges)
length(exclusive_edges_d12_d03)

exclusive_edges_d03_d06 <- setdiff(intersect(edges_d03, edges_d06), common_edges)
length(exclusive_edges_d03_d06)

#I use setdiff twice to find all the unqiue egdes in one networks and then print its value
exclusive_edges_d12 <- setdiff(setdiff(edges_d12, edges_d06), edges_d03)
length(exclusive_edges_d12)
exclusive_edges_d06 <- setdiff(setdiff(edges_d06, edges_d12), edges_d03)
length(exclusive_edges_d06)
exclusive_edges_d03 <- setdiff(setdiff(edges_d03, edges_d06), edges_d12)
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

