setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/scripts_postprocess')
source("my_functions.R")

library(ggplot2)
library(cowplot)
library(Metrics)
library(MLmetrics)
library(gridExtra)

## Récupère les fichiers résultats
# no photomorph caribu passive
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/nophotomorpho_1tige_leg_sky5_1t/')
nophoto_legume_active_1 <- read.table('toto_17111_l-egume.csv', sep=';',stringsAsFactors = FALSE)
nophoto_legume_active_2 <- read.table('toto_17112_l-egume.csv', sep=';',stringsAsFactors = FALSE)
nophoto_caribu_passive_1 <- read.table('outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
nophoto_caribu_passive_2 <- read.table('outputs_ratp_passive_17112_l-egume.csv', sep=',', stringsAsFactors = FALSE)


# no photomorph caribu active
setwd('C:/Users/mwoussen/cdd/mesocentre/legume_caribu/legume_nophoto_caribu_active')
nophoto_caribu_active_1 <- read.table('toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
nophoto_caribu_active_2 <- read.table('toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)

# photomorph caribu passive
setwd('C:/Users/mwoussen/cdd/mesocentre/legume_caribu/legume_photo_caribu_passive')
photo_legume_active_1 <- read.table('toto_17111_l-egume.csv', sep=';',stringsAsFactors = FALSE)
photo_legume_active_2 <- read.table('toto_17112_l-egume.csv', sep=';',stringsAsFactors = FALSE)
photo_caribu_passive_1 <- read.table('outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
photo_caribu_passive_2 <- read.table('outputs_ratp_passive_17112_l-egume.csv', sep=',', stringsAsFactors = FALSE)

# photomorph caribu active
setwd('C:/Users/mwoussen/cdd/mesocentre/legume_caribu/legume_photo_caribu_active')
photo_caribu_active_1 <- read.table('toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_caribu_active_2 <- read.table('toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)


#### EPSI COUVERT
# légende des modèles
data_names <- list("default", "caribu passif", "caribu actif")

# situation no photomorpho ,no ramification, 1 shoot
dataframes_1 <- list(nophoto_legume_active_1, nophoto_caribu_passive_1, nophoto_caribu_active_1)
dataframes_2 <- list(nophoto_legume_active_2, nophoto_caribu_passive_2, nophoto_caribu_active_2)

# situation photomorpho ,ramification, 3 shoots
dataframes_1 <- list(photo_legume_active_1, photo_caribu_passive_1, photo_caribu_active_1)
dataframes_2 <- list(photo_legume_active_2, photo_caribu_passive_2, photo_caribu_active_2)

# graphe
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
  xlab("Jour de l'année") +
  ylab(ytitle) +
  labs(color="Modèle")
p <- p+theme(plot.title = element_text(size=11))
p

#### LAI COUVERT
# légende des modèles
data_names <- list("default", "caribu actif")

# situation no photomorpho ,no ramification, 1 shoot
dataframes_1 <- list(nophoto_legume_active_1, nophoto_caribu_active_1)
dataframes_2 <- list(nophoto_legume_active_2, nophoto_caribu_active_2)

# situation photomorpho ,ramification, 3 shoots
dataframes_1 <- list(photo_legume_active_1, photo_caribu_active_1)
dataframes_2 <- list(photo_legume_active_2, photo_caribu_active_2)

# graphe
varname <- "SurfPlante"
ytitle <- "LAI"
globaltitle <- "LAI sur tout le couvert"
start <- 60
end <- 180
surfsol <- 0.16
df_1 <- plot_dataframe_canopy(varname, dataframes_1, data_names)
df_2 <- plot_dataframe_canopy(varname, dataframes_2, data_names)
df <- df_1
df$values <- df$values + df_2$values
df$values <- sapply(df$values, function(x){x/surfsol})
p <- ggplot(data=df, aes(x=x, y=values, group=light)) +
  geom_line(aes(color=light)) +
  ggtitle(globaltitle) +
  xlab("Jour de l'année") +
  ylab(ytitle) +
  labs(color="Modèle")
p <- p+theme(plot.title = element_text(size=11))
p

#### CARIBU%L-EGUME sur le epsi
# situation no photomorpho ,no ramification, 1 shoot
dataframes_1 <- list(nophoto_legume_active_1, nophoto_caribu_passive_1)
dataframes_2 <- list(nophoto_legume_active_2, nophoto_caribu_passive_2)

# situation photomorpho ,ramification, 3 shoots
dataframes_1 <- list(photo_legume_active_1, photo_caribu_active_1)
dataframes_2 <- list(photo_legume_active_2, photo_caribu_active_2)

df <- df_correlation_plante("PARiPlante", list(nophoto_legume_active_1,
                                         nophoto_legume_active_2,
                                         photo_legume_active_1,
                                         photo_legume_active_2), list(nophoto_caribu_passive_1,
                                                                      nophoto_caribu_passive_2,
                                                                      photo_caribu_passive_1,
                                                                      photo_caribu_passive_2), list("photomorpho OFF",
                                                                                                     "photomorpho OFF",
                                                                                                     "photomorpho ON",
                                                                                                     "photomorpho ON"))
#df <- df[rowSums(df[,1:2]) >0,]

# df <- df[df$lightmodel < 10,]

df_xy <- data.frame(x=c(0,0.16),y=c(0,0.16))

p <- ggplot() +
  geom_point(data=df, aes(x=default, y=lightmodel, group=situation, color=situation, shape=situation),alpha=0.5) +
  geom_line(data=df_xy,aes(x=x,y=y)) +
  labs(x="Défaut l-egume",y="CARIBU Passif",title="PAR intercepté par plante en W.m-2")
p

rmse_on <- rmse(df[df$situation=="photomorpho ON",]$default, df[df$situation=="photomorpho ON",]$lightmodel)
rmse_off <- rmse(df[df$situation=="photomorpho OFF",]$default, df[df$situation=="photomorpho OFF",]$lightmodel)

mae_off <- mae(df[df$situation=="photomorpho OFF",]$default, df[df$situation=="photomorpho OFF",]$lightmodel)
mae_on <- mae(df[df$situation=="photomorpho ON",]$default, df[df$situation=="photomorpho ON",]$lightmodel)


df_rmse <- data.frame(Situation=c("photomorpho ON", "photomorpho OFF"), RMSE=c(rmse_on, rmse_off), MAE=c(mae_on, mae_off))

t <- tableGrob(df_rmse, rows = NULL)

plot_grid(p, t, nrow=2, rel_heights = c(1, 0.5))


##### LAI PAR PLANTE
# library(RPackVGL)
#
# # préparation de la dataframe finale
# ls_toto <- c('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_caribu/nophotomorph_noramif_1shoots_cari_act/toto_17111_l-egume.csv',
#              'C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_caribu/nophotomorph_noramif_1shoots_cari_act/toto_17112_l-egume.csv')
# my_ltoto <- read_ltoto(ls_toto)
# names(my_ltoto)
#
#
# # CARIBU ACTIVE sans photomorpho
# simmoy <- build_simmoy(my_ltoto, lsusm=names(my_ltoto))
# simsd <- build_simmoy(my_ltoto, lsusm=names(my_ltoto), optSD = T)
#
# gg_plotsim("LAI", simmoy, simsd, name="Par plante, sans photomorphogénèse")

####### moyenne des LAI des plantes pour chaque pas de temps
dataframe_list = list(nophoto_caribu_active_1, nophoto_caribu_active_2)
values <- list()
steps <- list()
variable = "SurfPlante"
surfsol <- 0.16
for (i in 1:length(dataframe_list))
{
  for (j in 1:32)
  {
    values <- append(values, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,j+2], function(x) as.numeric(as.character(x))))
    steps <- append(steps, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,2], function(x) as.numeric(as.character(x))))
  }
}
df <- data.frame(x=sapply(steps,c), values=sapply(values,c))

values <- list()
steps <- list()
situation <- list()
for (i in 60:178)
{
  values <- append(values, sapply(mean(df[df$x == i ,2])/surfsol, function(x) as.numeric(as.character(x))))
  steps <- append(steps, as.numeric(i))
  situation <- append(situation, "caribu")
}

df_moy <- data.frame(x=sapply(steps,c), values=sapply(values,c), model=sapply(situation, c))

dataframe_list = list(nophoto_legume_active_1, nophoto_legume_active_2)
values <- list()
steps <- list()
for (i in 1:length(dataframe_list))
{
  for (j in 1:32)
  {
    values <- append(values, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,j+2], function(x) as.numeric(as.character(x))))
    steps <- append(steps, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,2], function(x) as.numeric(as.character(x))))
  }
}
df <- data.frame(x=sapply(steps,c), values=sapply(values,c))

values <- list()
steps <- list()
situation <- list()
for (i in 60:178)
{
  values <- append(values, sapply(mean(df[df$x == i ,2])/surfsol, function(x) as.numeric(as.character(x))))
  steps <- append(steps, as.numeric(i))
  situation <- append(situation, "legume")
}

df_moy_leg <- data.frame(x=sapply(steps,c), values=sapply(values,c), model=sapply(situation, c))

df_moy <- rbind(df_moy, df_moy_leg)

###### sd
dataframe_list = list(nophoto_caribu_active_1, nophoto_caribu_active_2)
values <- list()
steps <- list()
for (i in 1:length(dataframe_list))
{
  for (j in 1:32)
  {
    values <- append(values, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,j+2], function(x) as.numeric(as.character(x))))
    steps <- append(steps, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,2], function(x) as.numeric(as.character(x))))
  }
}
df <- data.frame(x=sapply(steps,c), values=sapply(values,c))

values <- list()
steps <- list()
situation <- list()
for (i in 60:178)
{
  values <- append(values, sapply(sd(df[df$x == i ,2])/surfsol, function(x) as.numeric(as.character(x))))
  steps <- append(steps, as.numeric(i))
  situation <- append(situation, "caribu")
}

df_sd <- data.frame(x=sapply(steps,c), sd=sapply(values,c), model=sapply(situation, c))

dataframe_list = list(nophoto_legume_active_1, nophoto_legume_active_2)
values <- list()
steps <- list()
for (i in 1:length(dataframe_list))
{
  for (j in 1:32)
  {
    values <- append(values, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,j+2], function(x) as.numeric(as.character(x))))
    steps <- append(steps, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,2], function(x) as.numeric(as.character(x))))
  }
}
df <- data.frame(x=sapply(steps,c), values=sapply(values,c))

values <- list()
steps <- list()
situation <- list()
for (i in 60:178)
{
  values <- append(values, sapply(sd(df[df$x == i ,2])/surfsol, function(x) as.numeric(as.character(x))))
  steps <- append(steps, as.numeric(i))
  situation <- append(situation, "legume")
}

df_sd_leg <- data.frame(x=sapply(steps,c), sd=sapply(values,c), model=sapply(situation, c))

df_sd <- rbind(df_sd, df_sd_leg)

df <- df_moy
df$sd <- df_sd$sd
df$ymin <- df$values - df$sd
df[df$ymin < 0,] <- 0

# ymin <- list()
# ymax <- list()
# steps <- list()
# situation <- list()
# for (i in 60:178)
# {
#   for (n in c("caribu", "legume"))
#   {
#     ymin <- append(ymin, sapply(df_moy[df_moy$x==i && df_moy$model==n,2] - df_sd[df_sd$x==i && df_sd$model==n,2], function(x) as.numeric(as.character(x))))
#     ymax <- append(ymax, sapply(df_moy[df_moy$x==i && df_moy$model==n,2] + df_sd[df_sd$x==i && df_sd$model==n,2], function(x) as.numeric(as.character(x))))
#     steps <- append(steps, as.numeric(i))
#     situation <- append(situation, n)
#   }
# }
#
# df_sd <- data.frame(x=sapply(steps,c), ymin=sapply(ymin,c), ymax=sapply(ymax,c), model=sapply(situation, c))

### PLOT
min <- 0
max <- 1.5*max(df_moy$values)

plot_var <- ggplot(data=df, aes(x=x, y=values, group=model)) +
  geom_ribbon(aes(ymin=values-sd,ymax=values+sd, fill=model), alpha=0.15) +
  geom_line(aes(color=model),linewidth=1)+
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=8))+
  labs(title = "Moyenne et écart-type du LAI sur les plantes, sans photomorphogénèse",x = "DOY", y = "LAI", color="Modèle", fill="Modèle")

plot_var


#### CPUTIME % STEPS
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/nophotomorpho_1tige_ratp_sky5_1t/')
nophoto_ratp_active_1 <- read.table('toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
nophoto_ratp_active_2 <- read.table('toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)

data_names <- list("l-egume default", "caribu actif", "ratp actif")

# situation no photomorpho ,no ramification, 1 shoot
dataframes_1 <- list(nophoto_legume_active_1, nophoto_caribu_active_1, nophoto_ratp_active_1)
dataframes_2 <- list(nophoto_legume_active_2, nophoto_caribu_active_2, nophoto_ratp_active_2)

varname <- "time"
ytitle <- "CPU time en s"
globaltitle <- "Temps de simulation totale"
start <- 60
end <- 180
df_1 <- plot_dataframe_canopy(varname, dataframes_1, data_names)
df_2 <- plot_dataframe_canopy(varname, dataframes_2, data_names)
df <- df_1
df$values <- df$values + df_2$values
df$values <- df$values/32
p <- ggplot(data=df, aes(x=x, y=values, group=light)) +
  geom_line(aes(color=light)) +
  ggtitle(globaltitle) +
  xlab("Jour de l'année") +
  ylab(ytitle) +
  labs(color="Modèle")
p <- p+theme(plot.title = element_text(size=11))
p

