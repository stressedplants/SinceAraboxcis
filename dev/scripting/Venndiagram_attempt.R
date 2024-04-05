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

edges_d12 = paste(d12Net[,1], d12Net[,2], sep = '_')
edges_d06 = paste(d06Net[,1], d06Net[,2], sep = '_')
edges_d03 = paste(d03Net[,1], d03Net[,2], sep = '_')

common_edges <- intersect(intersect(edges_d12, edges_d06), edges_d03)
length(common_edges)

exclusive_edges_d12_d06 <- setdiff(intersect(edges_d12, edges_d06), common_edges)
length(exclusive_edges_d12_d06)
exclusive_edges_d12_d03 <- setdiff(intersect(edges_d12, edges_d03), common_edges)
length(exclusive_edges_d12_d03)
exclusive_edges_d03_d06 <- setdiff(intersect(edges_d03, edges_d06), common_edges)
length(exclusive_edges_d03_d06)

exclusive_edges_d12 <- setdiff(setdiff(edges_d12, edges_d06), edges_d03)
length(exclusive_edges_d12)
exclusive_edges_d06 <- setdiff(setdiff(edges_d06, edges_d12), edges_d03)
length(exclusive_edges_d06)
exclusive_edges_d03 <- setdiff(setdiff(edges_d03, edges_d06), edges_d12)
length(exclusive_edges_d03)


d03_data <- c(length(exclusive_edges_d03), length(exclusive_edges_d03_d06), length(exclusive_edges_d12_d03), length(common_edges))
d06_data <- c(length(exclusive_edges_d03_d06), length(exclusive_edges_d06), length(exclusive_edges_d12_d06), length(common_edges))
d12_data <- c(length(exclusive_edges_d12_d03), length(exclusive_edges_d12_d06), length(exclusive_edges_d12), length(common_edges))

Venn_matrix <- matrix(c(d03_data, d06_data, d12_data), nrow = 3, ncol = 4, byrow = TRUE, 
                           dimnames = list(c("d03", "d06", "d12"), 
                                           c("shared_d03", "shared_d06", "shared_d12", "common")))

Venn_matrix <- cbind(Venn_matrix, rowSums(Venn_matrix))

# Update column names
colnames(Venn_matrix) <- c(colnames(Venn_matrix)[1:4], "total")

# Print the updated matrix
print(Venn_matrix)

