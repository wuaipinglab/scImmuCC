
# load("./data/test_data.rda",package="scImmuCC")
data("./data/test_data.rda",package="scImmuCC")
test_data <- as.matrix(test_data)
test <- scImmuCC_Layered(test_data,Non_Immune=FALSE)
