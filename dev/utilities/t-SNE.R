#tSNE
#had a go at doing t-SNE visulatiation, just change the code to fit your data 

#install all the required packages 
install.packages("Rtsne")
install.packages("ggplot2")

#load the required packages 
library(Rtsne)
library(ggplot2)

#had to covert gbox_filterd into a matrix format using as.matrix()
gbox_matrix <- as.matrix(gbox_filtered)

#Use the Rtsne function to perform t-SNE on your data. Specify the appropriate 
#parameters such as perplexity and dims (number of dimensions) based on your data:
tsne_result <- Rtsne(gbox_matrix, dims = 2, perplexity = 45)

#plot t-SNE
plot(tsne_result$Y, col =colours[clust[includeCells]], pch = 20, main = "t-SNE Plot Seedling_3d")
