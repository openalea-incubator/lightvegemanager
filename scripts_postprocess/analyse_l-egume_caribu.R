setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/scripts_postprocess')
source("my_functions.R")

library(ggplot2)
library(cowplot)
library(Metrics)
library(MLmetrics)
library(gridExtra)

## R?cup?re les fichiers r?sultats
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_caribu/nophotomorph_noramif_1shoots_cari_pass')
legume_active_1 <- read.table('toto_17111_l-egume.csv', sep=';',stringsAsFactors = FALSE)
legume_active_2 <- read.table('toto_17112_l-egume.csv', sep=';',stringsAsFactors = FALSE)
caribu_passive_1 <- read.table('outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
caribu_passive_2 <- read.table('outputs_ratp_passive_17112_l-egume.csv', sep=',', stringsAsFactors = FALSE)
dataframes_1 <- list(legume_active_1, caribu_passive_1)
dataframes_2 <- list(legume_active_2, caribu_passive_2)
data_names <- list("default", "caribu passive")


#### epsi
varname <- "epsi"
ytitle <- "Epsi"
globaltitle <- "Somme du epsilon sur tout le couvert"
start <- 60
end <- 180
df_1 <- plot_dataframe_canopy(varname, dataframes_1, data_names)
df_2 <- plot_dataframe_canopy(varname, dataframes_2, data_names)
df <- df_1
df$values <- df$values + df_2$values
p <- ggplot(data=df, aes(x=x, y=values, group=light)) +
  geom_line(aes(color=light)) +
  ggtitle(globaltitle) +
  xlab("Jour de l'ann?e") +
  ylab(ytitle) +
  labs(color="Mod?le")
p <- p+theme(plot.title = element_text(size=11))
p

#### LAI
varname <- "SurfPlante"
ytitle <- "LAI"
globaltitle <- "LAI sur tout le couvert"
start <- 60
end <- 180
surfsol <- 0.16
df <- plot_dataframe_canopy(varname, dataframes, data_names)
df$values <- sapply(df$values, function(x){x/surfsol})
p <- ggplot(data=df, aes(x=x, y=values, group=light)) +
  geom_line(aes(color=light)) +
  ggtitle(globaltitle) +
  xlab("Jour de l'ann?e") +
  ylab(ytitle) +
  labs(color="Mod?le")
p <- p+theme(plot.title = element_text(size=11))
p

