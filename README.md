# sc-ImmuCC: Hierarchical annotation for immune cell types in single-cell RNA-Seq 

The below demonstrates the result of Hierarchical annotation for some immune cells in the E-MTAB-11536 dataset (included in the package)
<div style="backgroud-color: #f5f5f5; padding: 10px">
library(scImmuCC)

data(package="scImmuCC)
data(test_data,package="scImmuCC") # load the test data

count <- as.matrix(test_data) # Convert test data to matrix

test <- scImmuCC_Layered(count = count ,Non_Immune = FALSE)

</div>

#In


