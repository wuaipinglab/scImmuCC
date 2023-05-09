#' @docType data
#' @name genelist
#' @format Data frames with cell type as column names.



layer0_genelist <- read.csv(file=paste("./extdata/","immune_non-immune_genelist(layer0).csv",sep=""),header=T)
layer0_genelist

layer1_genelist <- read.csv(file=paste("./extdata/","main_type_layer1_genelist.csv",sep=""),header=T)
layer1_genelist


Tcell_genelist <- read.csv(file=paste("./extdata/","Tcell_layer2_genelist.csv",sep=""),header=T)
Tcell_genelist

Bcell_genelist <- read.csv(file=paste("./extdata/","Bcell_layer2_genelist.csv",sep=""),header=T)
Bcell_genelist

DC_genelist <- read.csv(file=paste("./extdata/","DC_layer2_genelist.csv",sep=""),header=T)
DC_genelist

NK_genelist <- read.csv(file=paste("./extdata/","NK_layer2_genelist.csv",sep=""),header=T)
NK_genelist

Monocyte_genelist <- read.csv(file=paste("./extdata/","Monocyte_layer2_genelist.csv",sep=""),header=T)
Monocyte_genelist

Macrophage_genelist <- read.csv(file=paste("./extdata/","Macrophage_layer2_genelist.csv",sep=""),header=T)
Macrophage_genelist

ILC_genelist <- read.csv(file=paste("./extdata/","ILC_layer2_genelist.csv",sep=""),header=T)
ILC_genelist

CD4_genelist <- read.csv(file=paste("./extdata/","CD4_layer3_genelist.csv",sep=""),header=T)
CD4_genelist

CD8_genelist <- read.csv(file=paste("./extdata/","CD8_layer3_genelist.csv",sep=""),header=T)
CD8_genelist


# usethis::use_data(layer1_genelist)
# usethis::use_data(layer0_genelist)
# usethis::use_data(Tcell_genelist)
# usethis::use_data(Bcell_genelist)
# usethis::use_data(DC_genelist)
# usethis::use_data(NK_genelist)
# usethis::use_data(Monocyte_genelist)
# usethis::use_data(Macrophage_genelist)
# usethis::use_data(ILC_genelist)
# usethis::use_data(CD4_genelist)
# usethis::use_data(CD8_genelist)

