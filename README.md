# sc-ImmuCC
Hierarchical annotation of immune cells in scRNA-Seq data based on ssGSEA algorithm.

The below demonstrates the result of Hierarchical annotation for some immune cells in the E-MTAB-11536 dataset (included in the package)

library(scImmuCC)

data(package="scImmuCC)
data(test_data,package="scImmuCC") # load the test data

count <- as.matrix(test_data) # Convert test data to matrix

test <- scImmuCC_Layered(count = count ,Non_Immune = FALSE)
