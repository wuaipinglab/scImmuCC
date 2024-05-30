#' @title scImmuCC_main.
#' @description Calculation of gene enrichment scores and annotate cell type .
#' @details Input takes a count, genelist and filename, returns a data frame with cell annotation .
#' @param count  a matrix with cell unique barcodes as column names and gene names as row names .
#' @param genematrix  a data frame with cell types as column names .
#' @param ssGSEA_result  a data frame , scImmuCC return result .
#' @param filename  custom file name， character .
#' @return a data frame with cell annotation.
#' @import GSVA
#' @importFrom GSVA gsva



scImmuCC_main <- function(count,genematrix,filename){


  genelist <- as.list(genematrix)
  genelist <- lapply(genelist,function(x) x[!is.na(x)])
                     
  ssgsea <- ssgseaParam(
    count,
    genelist
    )
    #,
    #assay = NA_character_,
    #annotation = NA_character_,
    #minSize = 1,
    #maxSize = Inf,
    #alpha = 0.25,
    #normalize = TRUE
  #)
  
  ssgsea_score <- gsva(ssgsea)
  ##ssgsea_score = gsva(count, genelist, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
  score <- t(ssgsea_score)

  barcodes <- c()
  celltype <- c()
  for(i in 1:length(score[,1])){
    a1 <- score[i,]
    b1 <- as.data.frame(a1)
    b1$celltype <- rownames(b1)
    f <- b1[which(b1[,1]==max(b1[,1])),]
    barcodes <- c(barcodes,rownames(score)[i])
    celltype <- c(celltype,rownames(f))
  }

  ssGSEA_result <- data.frame(barcodes=barcodes,cell_type=celltype)
  head(ssGSEA_result)
  write.csv(ssGSEA_result,paste(filename,"_scImmuCC_label.csv",sep=""))
  return(ssGSEA_result)
}
