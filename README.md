# sc-ImmuCC: Hierarchical annotation for immune cell types in single-cell RNA-Seq 


# Installation
R programming language >= 4.1.1 ，packages Seurat and GSVA are required to use scImmuCC.

The installation from GitHub is in experimental stage but it can be used normally when the dependencies are installed:
<div style="backgroud-color: #f5f5f5; padding: 10px">
    
    if (!requireNamespace("Seurat", quietly = TRUE))
        install.packages("Seurat")
    
    if (!requireNamespace("GSVA", quietly = TRUE))
        BiocManager::install("GSVA")

    if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")
    
    remotes::install_github("wuaipinglab/scImmuCC")
</div>

# Example
Below is an example of Hierarchical annotation for some immune cells in the E-MTAB-11536 dataset (included in the package)

<div style="backgroud-color: #f5f5f5; padding: 10px">
    
    library(scImmuCC)

    #data(package="scImmuCC)
    data(test_data,package="scImmuCC") # load the test data

    count <- as.matrix(test_data) # Convert test data to matrix

    test <- scImmuCC_Layered(count = count ,Non_Immune = FALSE)

</div>


# QuickStart
The following is a quick tutorial on how to use scImmuCC to annotate immune cell types in scRNA-Seq.

<div style="backgroud-color: #f5f5f5; padding: 10px">
    
    library(scImmuCC)
    
    count <- read.csv(file=filename)     ##read your scRNA-Seq file
             
    count <- as.matrix(count)     ##  a matrix with cell unique barcodes as column names and gene names as row names
    
    test <- scImmuCC_Layered(test_data,Non_Immune=FALSE)    ##if your data have non-immune cell, Non_Immune = TRUE
            
</div>
The annotation results will be output in your current running directory。
            
            










