library(Matrix)
source('dev/utilities/dataprocessingHelperFunctions.R')
#Load the single cell data for seed 12
data = load('data/GSE226097_seedling_12d_230221.RData')
#loading in the og AraBoxcis network 
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)


# Filtering gbox ------------------------------------------------


dim(gbox)
#remove genes that have less than % of the genes expressed at all. 
thresh = 0.01
numberGenesPerCell = apply(gbox, 2, function(i){length(which(i > 0))})
includeCells = which(numberGenesPerCell > (thresh*dim(gbox)[1]))

gbox_filtered = gbox[,includeCells]

dim(gbox_filtered)


# UMAP data generation and visulisation -------------------------

library(umap)

pca = prcomp(gbox_filtered, scale. = T, rank. = 5)
gbox.pca.umap <- umap(pca$x)

colours = rainbow(length(unique(clust)))

plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2],
     col = "blue", pch = 20, cex = 0.5,
     main = 'PCA UMAP day 12 Seedling', xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')


# Using monocle3 ------------------------------------------------

library(Seurat)
library(progress)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)
library(patchwork)
library(dplyr)
library(Signac)
library(Matrix)
library(irlba)

#making gbox_filtered into a matrix so it can be convered to a 
#monocle cell data set
gbox_matrix_filtered <- as.matrix(gbox_filtered)
cds <- new_cell_data_set(gbox_matrix_filtered)
cds <- preprocess_cds(cds, num_dim = 5)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)

cds <- cluster_cells(cds, reduction_method = ("UMAP"), k = 15)
cds <- learn_graph(cds)
plot_cells(cds)
plot_cells(cds, color_cells_by = "partition", label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5,
           cell_size = 0.2,)
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           group_cells_by = "partition",
           show_trajectory_graph = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5,
           cell_size = 0.2) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

num_cells <- nrow(cds)
print(num_cells)


UMAP_data <- cds@int_colData@listData$reducedDims$UMAP
UMAP_data <- as.data.frame(UMAP_data)
frequencyplot_clustering_seedling_d12 <- plot(table(sort(clust)),
                                              xlab = 'Cluster Name',
                                              ylab = 'Number of Cells',
                                              main = 'Seedling - Day 12')
colours = rainbow(length(unique(clust)))
valid_clusters <- clust[rownames(UMAP_data)]
UMAP_data$cluster <- valid_clusters

UMAP_data$cluster <- as.numeric(as.character(UMAP_data$cluster))


cluster_colors <- c("#203E4E", "#2E6069", "#4E9D9E", "#53A798", "#64C39F", "#80BF9D",
                    "#84BF97", "#BAD5A9", "#E6D496", "#E7BD65", "#E6AA41", "#DA8835",
                    "#C77B3A", "#C7702E", "#B3724A", "#C4827F", "#C584A2", "#E26581",
                    "#E2607D")

cluster_colors_heatmap <- c("grey", "yellow", "grey", "red", "grey", "grey",
                            "grey", "grey", "grey", "yellow", "yellow", "grey",
                            "grey", "yellow", "red", "blue", "blue", "blue",
                            "red")

cluster_colors_heatmap_minor <- c("grey", "grey", "grey", "grey", "grey", "grey",
                            "grey", "grey", "grey", "grey", "grey", "grey",
                            "grey", "grey", "grey", "#C4827F", "#C584A2", "#E26581",
                            "grey")

cluster_colors_grey <- c("grey", "grey", "grey", "grey", "grey", "grey",
                         "grey", "grey", "grey", "grey", "grey", "grey",
                         "grey", "grey", "grey", "#C4827F", "#C584A2", "#E26581",
                         "grey")


# Plot with reordered legend and orginal ATLAS colors
ggplot(UMAP_data, aes(x = V1, y = V2, color = reorder(cluster, cluster))) +
  geom_point(alpha = 1, size = 0.01) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  ggtitle("UMAP by Atlas clusters") +
  scale_color_manual(values = cluster_colors, guide = guide_legend(override.aes = list(size = 5))) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot with reordered legend and greyed out colors
ggplot(UMAP_data, aes(x = V1, y = V2, color = reorder(cluster, cluster))) +
  geom_point(alpha = 1, size = 0.01) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  ggtitle("UMAP by Atlas clusters") +
  scale_color_manual(values = cluster_colors_grey, guide = guide_legend(override.aes = list(size = 5))) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot with reordered legend and greyed out colors with heatmap trio clusters
ggplot(UMAP_data, aes(x = V1, y = V2, color = reorder(cluster, cluster))) +
  geom_point(alpha = 1, size = 0.01, show.legend = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  ggtitle("UMAP by Heatmap Clusters") +
  scale_color_manual(values = cluster_colors_heatmap) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")


pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))


# Finding genes that change as a function of pseudotime ---------

deg_cds <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

num_significant <- sum(deg_cds$q_value < 0.01)

pseudotime_genes <- deg_cds %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head(num_significant)

top_pseudotime_genes <- pseudotime_genes %>% 
  arrange(desc(morans_test_statistic)) %>% 
  head(50)

top_pseudotime_genes_names <- rownames(top_pseudotime_genes)

PT_stress_genes <- c("AT3G05890", "AT1G35720", "AT5G46870", "AT2G15970",
                   "AT4G11600", "AT3G05500", "AT4G40010", "AT1G67090",
                   "AT3G05020", "AT1G27330", "AT4G39090", "AT4G19230" )
Pt_Photo_cellwall_genes <- c("AT1G67090", "AT1G79040", "AT2G44620",
                             "AT3G02230", "AT5G15650")
figure_genes <- c("AT3G05890", "AT4G19230", "AT3G05020",
                  "AT5G15650", "AT5G40340", "AT1G79040")



# plotting gene expression --------------------------------------

library(ggplot2)

#loop attempt 
plot_list <- list()

# Iterate over each term in the vector PT_corr_genes
for (x in top_pseudotime_genes_names) {
  # Extract color values for the current term x
  color_values <- gbox_filtered[x, ]
  colors <- colorRampPalette(c("grey90", "red"))(length(unique(color_values)))
  
  # Add color information to UMAP_data
  UMAP_data$color_X <- colors[match(color_values, unique(color_values))]
  
  # Create the plot for the current term x
  plot <- ggplot(UMAP_data, aes(x = V1, y = V2, color = color_X)) +
    geom_point(alpha = 1, size = 0.01) +
    scale_color_identity() +
    ggtitle(x) +
    theme_void() 
  
  # Add the plot to the list
  plot_list[[x]] <- plot
}

library(patchwork)

# Combine all plots in plot_list into one big figure
combined_plot <- wrap_plots(plotlist = plot_list)

# Print the combined plot
print(combined_plot)


# Heat maps ------------------
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
class(simpleNames)

#renaming varible to id doublicates 
simple_name_pretty = c("Unannotated 1", "Unannotated 2", "Phloem", "Unannotated 3",
                       "Trichoblast, Trichome 1", "Unannotated 4",
                       "Unannotated 5", "Unannotated 6", "Guard 1", "Trichome 1",
                       "Trichome 2", "Trichoblast, Trichome 2", "Phloem",
                       "Mesophyll", "Guard 2", "Meristem (Root and Shoot)",
                       "Procambium_PP 1", "Procambium_PP 2", "Vascular" )

# Print the updated vector
print(simple_name_pretty)
length(simple_name_pretty)
length(simpleNames)

rownames(geneExpByCluster) = simple_name_pretty

pheatmap(geneExpByCluster, scale = 'column', 
         main = "Heatmap by cluster", 
         show_colnames = FALSE)

#saving the underlying table as a table. 
write.table(geneExpByCluster, 
            'data/Seedling_d12_avgExpressionByCluster.txt')

#Focusing on transcription factor 
tfs = unique(araboxcis[,1])
tfSubs = tfs[which(tfs %in% rownames(gbox))]
pheatmap(geneExpByCluster[,tfSubs], scale = 'column', 
         main = "Heatmap by clustered cell type, Seedlings d12")

#finding correlation between every TF and every potential target
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
net_5 = GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees = 5) 

save(net, file = 'seedling_d12_network_nTree_5.RData')

ginieOutput = convertToAdjacency(net_50, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]


# Overlap -------------------------------------------------------



a = load("data/seedling_d12_network_nTree_50.RData")
newNet = GENIE3::getLinkList(net_50) #converts the weight matrix 
#(the gene network) into a sorted list of regulatory links. 
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)

#finding the size of my network
#sort by weight
newNet <- newNet[order(-newNet$weight), ]
#plot
plot(newNet$weight, main = "Weight vs Row Index", 
     xlab = "Row Index", ylab = "Weight",
     pch = 20,  # Smaller points, using the solid circle symbol
     col = "blue",
     cex = 0.6)  # Smaller point size

# Customize the y-axis with more divisions
axis(side = 2, at = pretty(newNet$weight, n = 10))

# Add a grid for better readability
grid()

if (!require(plotly)) {
  install.packages("plotly")
}

# Load plotly library
library(plotly)

# Create the interactive plot
interactive_plot <- plot_ly(data = newNet, 
                            x = ~seq_len(nrow(newNet)),  # Row index for x-axis
                            y = ~weight,  # Weight for y-axis
                            type = 'scatter',
                            mode = 'markers',
                            marker = list(size = 6, color = 'blue')) %>%
  layout(title = "Weight vs Row Index",
         xaxis = list(title = "Row Index"),
         yaxis = list(title = "Weight"))

# Display the interactive plot
interactive_plot

#filter network for weight with less than 0.014 - based on graph
newNet_filtered <- newNet[newNet$weight >= 0.014, ]
#size = 24155

#filter network for weight with less than 0.014 - based on araboxcis lowest weight
newNet_filtered_ara <- newNet[newNet$weight >= 0.008714, ]
#size = 81665

#finding overlaps between the single cell gene netork and AraBOXcis network

#gets the set of unique genes in your new network
genesInNet = unique(c(newNet[,1], newNet[,2]))

length(genesInNet)

#filter the AraBOXcis network to only contain genes that are in your new network
araboxcisFiltered = araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]

#extract the top edges in your new network, to make your network the same size as the araboxcisFiltered network.
newNetTopEdges = newNet[1:length(araboxcisFiltered[,1]),]

#reformat edges so it is more straightforward to compare them
edgesNew = paste(newNetTopEdges[,1], newNetTopEdges[,2], sep = '_')
edgesOld = paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep = '_')

#Now, you can come up with all the different parts of a Venn Diagram:

#The overlap

d12_araboxcis_overlap <- intersect(edgesNew, edgesOld)
length(which(edgesNew %in% edgesOld))
length(d12_araboxcis_overlap)

#in new network only
length(which(!(edgesNew %in% edgesOld)))

#in araBOXcis only
length(which(!(edgesOld %in% edgesNew)))


d12_araboxcis_overlap_genes <- strsplit(d12_araboxcis_overlap , "_")
regulatoryGene <- sapply(d12_araboxcis_overlap_genes, "[[", 1)
targetGene <- sapply(d12_araboxcis_overlap_genes, "[[", 2)
d12_araboxcis_overlap <- data.frame(regulatoryGene, targetGene)
d12_araboxcis_overlap <- merge(d12_araboxcis_overlap, newNet,
                               by.x = c("regulatoryGene", "targetGene"),
                               by.y = c("regulatoryGene", "targetGene"),
                               all.x = TRUE)

d12_araboxcis_overlap <- d12_araboxcis_overlap[order(-d12_araboxcis_overlap$weight), ]
d12_araboxcis_overlap_Net <- d12_araboxcis_overlap[1:1000, ]
write.csv(d12_araboxcis_overlap_Net,
          "data/d12_araboxcis_overlap_Net.cvs", row.names = FALSE)

d12_Net <- newNet[1:1000, ]
write.csv(d12_araboxcis_overlap_Net,
          "data/d12_Net.cvs", row.names = FALSE)

#comparing d12_Net and d12_araboxcis_overlap_net
edges_d12 = paste(d12_Net[,1], d12_Net[,2], sep = '_')
edges_d12_ara = paste(d12_araboxcis_overlap_Net[,1], d12_araboxcis_overlap_Net[,2], sep = '_')

net_edges_d12_araboxcis_overlap <- intersect(edges_d12, edges_d12_ara)
length(net_edges_d12_araboxcis_overlap)

# Finding important genes ---------------------------------------

tfsNew = table(newNetTopEdges[,1])
tfsOld = table(araboxcisFiltered[,1])[names(tfsNew)]

#histogram of degrees should look like an exponential distribution, 
#because biological networks are a kind of network called 'scale-free', 
#meaning that most TFs only regulate a small number of genes, but a 
#few are influential hubs.

hist(as.numeric(tfsNew), main = 'SinceAraBOXcis',
     xlab = 'degree of TFs')

hist(as.numeric(tfsOld), main = 'AraBOXcis',
     xlab = 'degree of TFs')

#histogram is not as expected (exponential distribution)



# Gene scores ---------------------------------------------------

library(igraph)
library(pheatmap)
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))
simple_network_ara <- graph_from_edgelist(as.matrix(araboxcisFiltered[,c(1,2)]))

#top 20 degrees genes
top_degrees <- sort(tfsNew, decreasing = TRUE)[1:20]
#top 20 degrees for araboxcis
top_degrees_ara <- sort(tfsOld, decreasing = TRUE)[1:20]

#top 20 betweeness genes
node_betweenness_all <- betweenness(simple_network)
node_betweenness = node_betweenness_all[which(node_betweenness_all > 0)]
top_betweenness <- sort(node_betweenness, decreasing = TRUE)[1:100]
node_betweenness <- as.data.frame(node_betweenness)

#top 20 betweeness genes for araboxcis
node_betweenness_all_ara <- betweenness(simple_network_ara)
node_betweenness_ara = node_betweenness_all_ara[which(node_betweenness_all_ara > 0)]
top_betweenness_ara <- sort(node_betweenness_ara, decreasing = TRUE)[1:100]
node_betweenness_ara <- as.data.frame(node_betweenness_ara)

#top 20 hub genes
node_hub_all <- hub_score(simple_network)$vector
node_hub = node_hub_all[which(node_hub_all > 0)]
top_hub <- sort(node_hub, decreasing = TRUE)[1:20]

top_betweenness <- as.data.frame(top_betweenness)
betweenness_genes <- rownames(top_betweenness)

top_betweenness_ara <- as.data.frame(top_betweenness_ara)
betweenness_genes_ara <- rownames(top_betweenness_ara)

top_hub <- as.data.frame(top_hub)
hub_genes <- rownames(top_hub)
degrees_top <- as.data.frame(top_degrees)
degrees_top_ara <- as.data.frame(top_degrees_ara)
degree_genes <- degrees_top[,1]
degree_genes_ara <- degrees_top_ara[,1]

#comparing degrees
top_genes <- data.frame(degree_genes, betweenness_genes, hub_genes)
colnames(top_genes) <- c("Degrees", "Betweenness", "Hub")

comparison_degree_genes <- data.frame(degree_genes, degree_genes_ara)
colnames(comparison_degree_genes) <- c("SinceAraBOXcis", "AraBOXcis")

common_strings <- intersect(degrees_top$Var1, degrees_top_ara$Var1)
print(common_strings)
#of top 20 degrees 5 are common.

#comaparing betweenness


  


install.packages("writexl")
library(writexl)

# Write the data frame to the Excel file
write_xlsx(top_genes, "data/top_genes.xlsx")

#graph of top degrees ara vs Scara

library(ggplot2)
library(ggrepel)

df_tfsNew <- as.data.frame(tfsNew)
colnames(df_tfsNew) <- c("genename", "degrees_tfsNew")
df_tfsOld <- as.data.frame(tfsOld)
colnames(df_tfsOld) <- c("genename", "degrees_tfsOld")
merged_tfs <- merge(df_tfsNew, df_tfsOld, by = "genename", all = TRUE)

merged_tfs$point_color <- ifelse(merged_tfs$genename %in% c("AT5G08130", "AT3G19290", "AT4G36730", "AT5G46760", "AT3G17100"), "highlight", "default")

# Create the plot of degrees
ggplot(data = merged_tfs, aes(x = degrees_tfsNew, y = degrees_tfsOld, color = point_color)) +
  geom_point() +
  geom_text_repel(aes(label = genename), size = 3, 
                  fontface = ifelse(merged_tfs$genename %in% c("AT5G08130", "AT3G19290", "AT4G36730", "AT5G46760", "AT3G17100"), "bold", "plain"),
                  color = "black", segment.color = "lightgray") +
  labs(x = "SinceAraBOXcis", y = "AraBOXcis", title = "Comparison of Degrees") +  # Corrected title
  theme_classic() +
  scale_x_continuous(limits = c(0, 2500)) +
  scale_color_manual(values = c("default" = "grey70", "highlight" = "red")) +
  theme(legend.position = "none",  # Remove the legend
        plot.title = element_text(hjust = 0.5))  # Center the title

# Create the plot of betweenness 
node_betweenness_ara$RowNames <- rownames(node_betweenness_ara)
node_betweenness$RowNames <- rownames(node_betweenness)

# Merge the data frames based on the new column 'RowNames'
merged_btwn <- merge(node_betweenness, node_betweenness_ara, by = "RowNames", all = TRUE)

# Create the plot of degrees
ggplot(data = merged_btwn, aes(x = node_betweenness, y = node_betweenness_ara)) +
  geom_point() +
  labs(x = "SinceAraBOXcis", y = "AraBOXcis", title = "Comparison of Betweenness Scores") +
  geom_text_repel(aes(label = RowNames), size = 3, color = "black", segment.color = "lightgray") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

# Finding GO terms  ---------------------------------------------

a = load('data/functionalData.RData')
source('dev/utilities/dataprocessingHelperFunctions.R')

pafwayOut = pafway(GOconcat, newNet_filtered, unique(goOfInterest))
rownames(pafwayOut) = colnames(pafwayOut)

#Filter to only include rows and columns with at least one significant factor:
atLeastOneSigRow = which(apply(pafwayOut, 1,
                               function(i){length(which(i < 0.05))}) > 0)


atLeastOneSigCol = which(apply(pafwayOut, 2,
                               function(i){length(which(i < 0.05))}) > 0)

pafwayInterestingOnly = pafwayOut[atLeastOneSigRow, atLeastOneSigCol]

#Here, the terms associated with the columns are upstream of the terms associated with the rows.  
#The values correspond to p-values, so smaller is more significant.
pheatmap(pafwayInterestingOnly)

my_colors <- colorRampPalette(c("blue", "white"))(100)

# Apply log10 transformation to your data
log_data <- log10(pafwayInterestingOnly)

# Plot the heatmap
pheatmap(log_data,
         color = my_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
