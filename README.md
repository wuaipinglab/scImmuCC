# sc-ImmuCC: Hierarchical annotation for immune cell types in single-cell RNA-Seq 

The below demonstrates the result of Hierarchical annotation for some immune cells in the E-MTAB-11536 dataset (included in the package)
<div style="backgroud-color: #f5f5f5; padding: 10px">
library(scImmuCC)

    data(package="scImmuCC)
    data(test_data,package="scImmuCC") # load the test data

    count <- as.matrix(test_data) # Convert test data to matrix

    test <- scImmuCC_Layered(count = count ,Non_Immune = FALSE)

</div>

# Installation
R programming language >= 4.1.1 is required to use scImmuCC.

The installation from GitHub is in experimental stage but gives the newest feature:

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("wuaipinglab/scImmuCC")

# QuickStart
The following is a quick tutorial on how to use scImmuCC to annotate immune cell types in scRNA-Seq.

library(scImmuCC)

count <- read.csv(file=filename) #read your scRNA-Seq file

count <- as.matrix(count)

test <- scImmuCC_Layered(test_data,Non_Immune=FALSE) #if your data have non-immune cell, Nn_Immune = TRUE

The annotation results will be output in your current running directory。










