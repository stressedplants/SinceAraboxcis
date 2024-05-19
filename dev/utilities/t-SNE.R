#tSNE
#had a go at doing t-SNE visualisation and made a function, just change the code to fit your data 

perform_tSNE_visualisation <- function(gbox_filtered, colours, clust, includeCells, dims = 2, perplexity = 45, seed = NULL) {
  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
    # Install required packages if not already installed
  if (!requireNamespace("Rtsne", quietly = TRUE))
    install.packages("Rtsne")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
  
  # Load required packages
  library(Rtsne)
  library(ggplot2)
  
  # Convert gbox_filtered to matrix format
  gbox_matrix <- as.matrix(gbox_filtered)
  
  # Use the Rtsne function to perform t-SNE
  tsne_result <- Rtsne(gbox_matrix, dims = dims, perplexity = perplexity)
  
  # Check if colours, clust, and includeCells are provided
  if (missing(colours) || missing(clust) || missing(includeCells)) {
    stop("Arguments 'colours', 'clust', and 'includeCells' must be provided.")
  }
  
  # Plot t-SNE
  plot(tsne_result$Y, col = colours[clust[includeCells]], pch = 20, main = "t-SNE Plot Seedling_3d")
}


# Example usage
perform_tSNE_visualisation(gbox_filtered, colours, clust, includeCells, seed = 123)



