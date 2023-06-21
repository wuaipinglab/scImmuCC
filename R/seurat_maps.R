#' @title Seurat Object.
#' @description Creating dynamic Seurat Object to add annotation information and draw images.
#' @details Input takes a count, genelist and scImmuCC result, returns tsne, umap, Dotplot and pheatmaps.
#' @param count  a matrix with cell unique barcodes as column names and gene names as row names .
#' @param genematrix  a data frame with cell types as column names .
#' @param ssGSEA_result  a data frame , scImmuCC1.1 return result .
#' @param filename  custom file name， character .
#' @return 4 pictures.
#' @import ggplot2
#' @import Seurat
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData PercentageFeatureSet VariableFeatures RunPCA RunUMAP RunTSNE DimPlot DotPlot DoHeatmap ProjectDim JackStraw ScoreJackStraw



seurat_Heatmap <- function(count,genematrix,ssGSEA_result,filename){

  count <- count[, !duplicated(colnames(count))]
  seurat.data <- CreateSeuratObject(counts = count, project = filename)#, min.cells = 3, min.features = 200)
  seurat.data
  seurat.data[["percent.mt"]] <- PercentageFeatureSet(seurat.data, pattern = "^MT-")
  HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes_total,rownames(seurat.data@assays$RNA))

  HB.genes <- rownames(seurat.data@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  seurat.data[["percent.HB"]]<-PercentageFeatureSet(seurat.data,features=HB.genes)

  seurat.data <- NormalizeData(seurat.data, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.data <- FindVariableFeatures(seurat.data, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat.data)
  seurat.data <- ScaleData(seurat.data, features = all.genes)

  #Perform linear dimensional reductionPerform linear dimensional reduction
  counts <- seurat.data[["RNA"]]@counts
  cells <- length(counts[2,])
  if(cells>50){
    seurat.data <- RunPCA(seurat.data, features = VariableFeatures(object = seurat.data))
  }else{
    seurat.data <- RunPCA(seurat.data,npcs = (cells-1), features = VariableFeatures(object = seurat.data))
  }
  
  head(seurat.data@reductions$pca@cell.embeddings)
  head(seurat.data@reductions$pca@feature.loadings)
  seurat.data <- ProjectDim(object = seurat.data)
  seurat.data <- JackStraw(seurat.data, num.replicate = 100)
  seurat.data <- ScoreJackStraw(seurat.data, dims = 1:20)

  seurat.data <- RunUMAP(seurat.data, dims = 1:10)
  n <- length(count[2,])
  #seurat.data <- seurat.data[ ,!duplicated(colnames(seurat.data))]
  
  #if(n>91){
  #  seurat.data <- RunTSNE(seurat.data, dims = 1:10)
  #}
  #else{
  #  seurat.data <- RunTSNE(seurat.data, dims = 1:10, perplexity = 1)
  #}
  seurat.data <- RunTSNE(seurat.data, dims = 1:10)
  head(seurat.data@reductions$tsne@cell.embeddings)

  #Add annotation information to Seurat object
  seurat.data@meta.data$cell_type_pred <- ssGSEA_result[,2]

  pdf(paste(filename, "_tSNE", ".pdf", sep=""), width=12, height=10)
  p1 <- DimPlot(seurat.data, reduction = "tsne", group.by = "cell_type_pred",label = TRUE, pt.size=1)
  plot(p1)
  dev.off()
  pdf(paste(filename, "_UMAP", ".pdf", sep=""), width=12, height=10)
  p2 <- DimPlot(seurat.data, reduction = "umap", group.by = "cell_type_pred",label = TRUE, pt.size=1)
  plot(p2)
  dev.off()

  genelist <- as.list(genematrix)
  genelist <- lapply(genelist,function(x) x[!is.na(x)])
  gene <- c()
  for(i in 1:length(genelist)){
    gene <- c(gene,genelist[[i]])
  }
  length(gene)
  gene <- unique(gene)
  Y <- intersect(gene,rownames(count))
  length(Y)
  features <- Y
  pdf(paste(filename,"_UMAP_DotPlot",".pdf",sep=""), width=12, height=10)
  p3 <- DotPlot(seurat.data, features = features,group.by = "cell_type_pred") + RotatedAxis()
  plot(p3)
  dev.off()

  pdf(paste(filename, "_DoHeatmap.pdf", sep=""), width=12, height=10)
  p4 <- DoHeatmap(seurat.data, features = Y,group.by = "cell_type_pred") + NoLegend()
  plot(p4)
  dev.off()

}

