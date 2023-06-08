#' @title scImmuCC_Layered.
#' @description Creating Hierarchical_annotation for scRNA-Seq data immune cell.
#' @details Input takes a cells-genes matrix with cell unique barcodes as column names and gene names as row names and returns the cells annotation.
#' @param count  a matrix with cell unique barcodes as column names and gene names as row names .
#' @param Non_Immune Whether non-immune cells are included in the matrix.
#' @return Data frames with the barcodes and cell types, and some maps.
#' @import GSVA
#' @importFrom GSVA gsva
#' @import Matrix
#' @export
#' @examples test_data
#' library(GSVA)
#' library(Seurat)

scImmuCC_Layered <- function(count,Non_Immune=TRUE){

  # if (!requireNamespace("GSVA", quietly = TRUE)) {
    BiocManager::install("GSVA")
  # }
  # if (!requireNamespace("Seurat", quietly = TRUE)){
    install.packages("Seurat")
  # }
  data("./data/layer1_genelist.rda",package="scImmuCC")
  data("./data/Tcell_genelist.rda",package="scImmuCC")
  data("./data/Bcell_genelist.rda",package="scImmuCC")
  data("./data/DC_genelist.rda",package="scImmuCC")
  data("./data/NK_genelist.rda",package="scImmuCC")
  data("./data/Monocyte_genelist.rda",package="scImmuCC")
  data("./data/Macrophage_genelist.rda",package="scImmuCC")
  data("./data/ILC_genelist.rda",package="scImmuCC")
  data("./data/CD4_genelist.rda",package="scImmuCC")
  data("./data/CD8_genelist.rda",package="scImmuCC")
  data("./data/layer0_genelist.rda",package="scImmuCC")


  if(Non_Immune==TRUE){

    layer0_result <- scImmuCC_main(count,layer0_genelist,"Layer0")
    immune <- layer0_result[which(layer0_result[,2]=="Immune"),]
    count_immune <- count[,immune[,1]]
    ssGSEA_result <- scImmuCC_main(count_immune,layer1_genelist,"Layer1")
    seurat_result <- seurat_Heatmap(count,layer1_genelist,ssGSEA_result,"Layer1")

  }else{

    ssGSEA_result <- scImmuCC_main(count,layer1_genelist,"Layer1")
    seurat_result <- seurat_Heatmap(count,layer1_genelist,ssGSEA_result,"Layer1")
  }

  cell_type <- unique(ssGSEA_result[,2])

  if("Tcell" %in% cell_type){
    sub_ssGSEA_Tcell <- ssGSEA_result[which(ssGSEA_result[,2]=="Tcell"),]
    sub_count_Tcell <- count[,sub_ssGSEA_Tcell[,1]]
    sub_count_Tcell <- as.matrix(sub_count_Tcell)
    ssGSEA_Tcell <- scImmuCC_main(sub_count_Tcell,Tcell_genelist,"Layer2_Tcell")
    seurat_Tcell <- seurat_Heatmap(sub_count_Tcell,Tcell_genelist,ssGSEA_Tcell,"Layer2_Tcell")

    cell_type2 <- unique(ssGSEA_Tcell[,2])

    if("CD4_T" %in% cell_type){
      sub_ssGSEA_CD4 <- ssGSEA_Tcell[which(ssGSEA_Tcell[,2]=="CD4_T"),]
      sub_count_CD4 <- count[,sub_ssGSEA_CD4[,1]]
      sub_count_CD4 <- as.matrix(sub_count_CD4)
      ssGSEA_CD4 <- scImmuCC_main(sub_count_CD4,CD4_genelist,"Layer3_CD4")
      cells <- length(sub_count_CD4[2,])
      if(cells>50){
        seurat_CD4 <- seurat_Heatmap(sub_count_CD4,CD4_genelist,ssGSEA_CD4,"Layer3_CD4")
      }

    }


    if("CD8_T" %in% cell_type){
      sub_ssGSEA_CD8 <- ssGSEA_Tcell[which(ssGSEA_Tcell[,2]=="CD8_T"),]
      sub_count_CD8 <- count[,sub_ssGSEA_CD8[,1]]
      sub_count_CD8 <- as.matrix(sub_count_CD8)
      ssGSEA_CD8 <- scImmuCC_main(sub_count_CD8,CD8_genelist,"Layer3_CD4")
      cells <- length(sub_count_CD8[2,])
      if(cells>50){
        seurat_CD8 <- seurat_Heatmap(sub_count_CD8,CD8_genelist,ssGSEA_CD8,"Layer3_CD8")
      }

    }
  }

  if("Bcell" %in% cell_type){
    sub_ssGSEA_Bcell <- ssGSEA_result[which(ssGSEA_result[,2]=="Bcell"),]
    sub_count_Bcell <- count[,sub_ssGSEA_Bcell[,1]]
    sub_count_Bcell <- as.matrix(sub_count_Bcell)
    ssGSEA_Bcell <- scImmuCC_main(sub_count_Bcell,Bcell_genelist,"Layer2_Bcell")
    cells <- length(sub_count_Bcell[2,])
    if(cells>50){
       seurat_Bcell <- seurat_Heatmap(sub_count_Bcell,Bcell_genelist,ssGSEA_Bcell,"Layer2_Bcell")
    }
  }

  if("DC" %in% cell_type){
    sub_ssGSEA_DC <- ssGSEA_result[which(ssGSEA_result[,2]=="DC"),]
    sub_count_DC <- count[,sub_ssGSEA_DC[,1]]
    sub_count_DC <- as.matrix(sub_count_DC)
    ssGSEA_DC <- scImmuCC_main(sub_count_DC,DC_genelist,"Layer2_DC")
    cells <- length(sub_count_DC[2,])
    if(cells>50){
      seurat_DC <- seurat_Heatmap(sub_count_DC,DC_genelist,ssGSEA_DC,"Layer2_DC")
    }

  }

  if("NK" %in% cell_type){
    sub_ssGSEA_NK <- ssGSEA_result[which(ssGSEA_result[,2]=="NK"),]
    sub_count_NK <- count[,sub_ssGSEA_NK[,1]]
    sub_count_NK <- as.matrix(sub_count_NK)
    ssGSEA_NK <- scImmuCC_main(sub_count_NK,NK_genelist,"Layer2_NK")
    cells <- length(sub_count_NK[2,])
    if(cells>50){
      seurat_NK <- seurat_Heatmap(sub_count_NK,NK_genelist,ssGSEA_NK,"Layer2_NK")
    }
  }

  if("Monocyte" %in% cell_type){
    sub_ssGSEA_Mono <- ssGSEA_result[which(ssGSEA_result[,2]=="Monocyte"),]
    sub_count_Mono <- count[,sub_ssGSEA_Mono[,1]]
    sub_count_Mono <- as.matrix(sub_count_Mono)
    ssGSEA_Mono <- scImmuCC_main(sub_count_Mono,Monocyte_genelist,"Layer2_Monocyte")
    cells <- length(sub_count_Mono[2,])
    if(cells>50){
      seurat_Mono <- seurat_Heatmap(sub_count_Mono,Monocyte_genelist,ssGSEA_Mono,"Layer2_Monocyte")
    }
  }

  if("Macrophage" %in% cell_type){
    sub_ssGSEA_Mac <- ssGSEA_result[which(ssGSEA_result[,2]=="Macrophage"),]
    sub_count_Mac <- count[,sub_ssGSEA_Mac[,1]]
    sub_count_Mac <- as.matrix(sub_count_Mac)
    ssGSEA_Mac <- scImmuCC_main(sub_count_Mac,Macrophage_genelist,"Layer2_Macrophage")
    cells <- length(sub_count_Mac[2,])
    if(cells>50){
      seurat_Mac <- seurat_Heatmap(sub_count_Mac,Macrophage_genelist,ssGSEA_Mac,"Layer2_Macrophage")
    }
  }

  if("ILC" %in% cell_type){
    sub_ssGSEA_ILC <- ssGSEA_result[which(ssGSEA_result[,2]=="ILC"),]
    sub_count_ILC <- count[,sub_ssGSEA_ILC[,1]]
    sub_count_ILC <- as.matrix(sub_count_ILC)
    ssGSEA_ILC <- scImmuCC_main(sub_count_ILC,ILC_genelist,"Layer2_ILC")
    cells <- length(sub_count_ILC[2,])
    if(cells>50){
      seurat_ILC <- seurat_Heatmap(sub_count_ILC,ILC_genelist,ssGSEA_ILC,"Layer2_ILC")
    }
  }
}
