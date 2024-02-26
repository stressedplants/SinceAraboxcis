# Making a function for creating a Seurat object: 

SeuratUMAP_function <- function(x){
  if(!require("Seurat", character.only = TRUE)) {
    remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE) 
    library("Seurat",  character.only = TRUE) 
  }
  my.SeuratObject <- CreateSeuratObject(counts = x, project = "Araboxcis UMAP")
  my.SeuratObject <- NormalizeData(my.SeuratObject)
  all_genes <- rownames(my.SeuratObject)
  my.SeuratObject <- ScaleData(my.SeuratObject, features = all_genes)
  my.SeuratObject <- FindVariableFeatures(my.SeuratObject, selection.method = "vst", nfeatures = 2000)
  my.SeuratObject <- RunPCA(my.SeuratObject, features = VariableFeatures(object = my.SeuratObject))
  my.SeuratObject <- FindNeighbors(my.SeuratObject, dims = 1:20)
  my.SeuratObject <- FindClusters(my.SeuratObject, resolution = 0.5)
  Seurat_tSNE_05Res <- DimPlot(my.SeuratObject, group.by = "RNA_snn_res.0.5", label = T)
  Idents(my.SeuratObject) <- "RNA_snn_res.0.5"
  my.SeuratObject <- RunUMAP(my.SeuratObject, dims = 1:20)
  Seurat_UMAP <- DimPlot(my.SeuratObject, reduction = 'umap') 
  return(Seurat_UMAP)
}

