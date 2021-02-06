# EXTERNAL SOURCE ---------------------------------------------------------------------------------------

source("CCApackage.R")

# PARAMETERS --------------------------------------------------------------------------------------------

  data_X <- CLUSTER1_IFN
data_Y <- CLUSTER1_TOLL

var_x <- c("IRF8",	"GBP1",	"BST2",	"USP18",	"IFITM1",	"IFITM3",	"IFI30",	"TRIM25",	"GBP5",	"IFNG",	"TRIM8")
var_y <- c("IRF7",	"RPS6KA1",	"CTSL",	"MAPKAPK3")


list_param <- list(
    empty = FALSE,
    heatmap_clust_method = "ward.D2", #pesquisar
    heatmap_color_ramp = usecol(pal = pal_unikn, n = 7),
    heatmap_border_color = grey(0.4),
    heatmap_fontsize = 10,
    plotCCA_canonical_variate = 1, #primeira e segunda variavel canonica
    plotCCA_corr_min = 0.5, 
    plotCCA_x_title = "Infectado - IFN genes",
    plotCCA_y_title = "Infectado - TLR genes",
    plotCCA_short_x_title = "IFN", #variavel canonica name
    plotCCA_short_y_title = "TLR", #same
    plotCorCCA_n_canonical_variates = 2, #reactoma
    plotCorCCA_corr_min = 0.7,
    plotCorCCA_seed = 10, #qualquer número, organização do graph
    plotCorCCA_node_size = 5,
    plotCorCCA_label_size = 3,
    plots_x_color = "deepskyblue3",
    plots_y_color = "gray1",
    plots_file_extension = "pdf", #pdf gera mais rápido
    plots_width = 13,
    plots_height = 7,
    data_log_trans = FALSE, #se não estiver normalizado
    data_log_base = 2, #aplicando normalização em log base2
    data_scale_trans = TRUE,
    data_sheet1 = 1,
    data_sheet2 = 2,
    data_X = data_X[, var_x],
    data_Y = data_Y[, var_y]
  )

dim(data_X)
dim(data_Y)
summary(data_X)
summary(data_Y)

dados <- cbind(data_X,data_Y)
dim(dados)
library(pheatmap)
install.packages("textshape")
library(textshape)

S <- cor(dados)[46:(45+14), 1:45]
dim(S)
clustmat <- cluster_matrix(S, dim = "both", method = "ward.D2")
pheatmap(clustmat, cluster_rows = FALSE, cluster_cols = FALSE)

head(dados)
heatmap(as.matrix(dados))
class(dados)
# CCA ---------------------------------------------------------------------------------------------------

runCCA(list_param, publish = FALSE) # run CCA using custom values for parameters

runCCA(list_param, publish = TRUE) # generate files for publication
