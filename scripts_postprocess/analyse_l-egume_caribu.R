setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/scripts_postprocess')
source("my_functions.R")

library(ggplot2)
library(cowplot)
library(Metrics)
library(MLmetrics)
library(gridExtra)

## Récupère les fichiers résultats
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_caribu')
legume_active <- read.table('toto_17112_l-egume.csv', sep=';',stringsAsFactors = FALSE)
caribu_passive <- read.table('outputs_ratp_passive_17112_l-egume.csv', sep=',', stringsAsFactors = FALSE)
dataframes <- list(legume_active, caribu_passive)
data_names <- list("default", "caribu")


#### epsi
varname <- "epsi"
ytitle <- "Epsi"
globaltitle <- "Somme du epsilon sur tout le couvert"
start <- 60
end <- 180
p_grid <- plot_dataframe_canopy(varname, dataframes, data_names)
p_grid

#### LAI
varname <- "SurfPlante"
ytitle <- "LAI"
globaltitle <- "LAI sur tout le couvert"
start <- 60
end <- 180
surfsol <- 0.16
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=FALSE, lai_surfsol = surfsol)
p_grid

globaltitle <- "LAI sur tout le couvert, ZOOM"
start <- 160
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=FALSE, lai_surfsol = surfsol)
p_grid