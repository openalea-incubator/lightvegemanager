setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/scripts_postprocess')
source("my_functions.R")

library(ggplot2)
library(cowplot)


########################
#                      #
#     A NETTOYER       #
#                      #
########################




## TEMPS CPU ##
file_list <- list(
                  
                  'legume_ratp_leg_sky5_1t_cputime_2ent.csv',
                  'legume_ratp_leg_sky5_25t_cputime_2ent.csv',
                  'legume_ratp_leg_sky46_cputime_2ent.csv',
                  'legume_ratp_leg_sky100_cputime_2ent.csv',
                  'legume_ratp_leg_sky5_1t_cputime_1ent.csv',
                  'legume_ratp_leg_sky5_25t_cputime_1ent.csv',
                  'legume_ratp_leg_sky46_cputime_1ent.csv',
                  'legume_ratp_leg_sky100_cputime_1ent.csv'
                  )
situation_list <- list("RATP: sky5 | 1 tube", "RATP: sky5 | 25 tubes", "RATP: sky46 | 25 tubes", "RATP: sky100 | 25 tubes")
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/cputimes/')

values <- list()
situation <- list()
nb_entity <- list()
for (i in 1:8)
{
  df <- read.table(file_list[[i]], sep=',',stringsAsFactors = FALSE)
  
  values <- append(values, as.numeric(df$V3[2]))
  if (i < 5)
  {
    situation <- append(situation, situation_list[[i]])
    nb_entity <- append(nb_entity, "2 entités") 
  }
  else
  {
    situation <- append(situation, situation_list[[i-4]])
    nb_entity <- append(nb_entity, "1 entité")
  }
}
# cas particulier de l-egume
df <- read.table(file_list[[1]], sep=',',stringsAsFactors = FALSE)
values <- append(values, as.numeric(df$V1[2]))
situation <- append(situation, "RiRi (l-egume)")
nb_entity <- append(nb_entity, "2 entités") 

df <- read.table(file_list[[5]], sep=',',stringsAsFactors = FALSE)
values <- append(values, as.numeric(df$V1[2]))
situation <- append(situation, "RiRi (l-egume)")
nb_entity <- append(nb_entity, "1 entité") 

data_cputime <- data.frame(x=sapply(situation,c), values=sapply(values,c), n_entity=sapply(nb_entity, c))

# total du cpu time
p1 <- ggplot(data=data_cputime, aes(x=reorder(x, values), y=values, group=n_entity)) + 
  geom_point(shape=15) +
  geom_line(aes(color=n_entity), size=1.2) +
  xlab("Situation") +
  ylab("CPU time en s")+
  labs(color="nombre d'entités")+
  theme(axis.text.x = element_text(angle = 90))

# courbe ratio RATP/l-egume par situation
file_list <- list(
  
  'legume_ratp_leg_sky5_1t_cputime_2ent.csv',
  'legume_ratp_leg_sky5_25t_cputime_2ent.csv',
  'legume_ratp_leg_sky46_cputime_2ent.csv',
  'legume_ratp_leg_sky5_1t_cputime_1ent.csv',
  'legume_ratp_leg_sky5_25t_cputime_1ent.csv',
  'legume_ratp_leg_sky46_cputime_1ent.csv'
)
situation_list <- list("RATP: sky5 | 1 tube", "RATP: sky5 | 25 tubes", "RATP: sky46 | 25 tubes")
values <- list()
situation <- list()
nb_entity <- list()
for (i in 1:6)
{
  df <- read.table(file_list[[i]], sep=',',stringsAsFactors = FALSE)
  
  values <- append(values, as.numeric(df$V3[2]) / as.numeric(df$V1[2]))
  if (i < 4)
  {
    situation <- append(situation, situation_list[[i]])
    nb_entity <- append(nb_entity, "2 entités") 
  }
  else
  {
    situation <- append(situation, situation_list[[i-3]])
    nb_entity <- append(nb_entity, "1 entité")
  }
}

data_ratio <- data.frame(x=sapply(situation,c), values=sapply(values,c), n_entity=sapply(nb_entity, c))

# graphe
p2 <- ggplot(data=data_ratio, aes(x=reorder(x, values), y=values, group=n_entity)) + 
  geom_point(shape=15) +
  geom_line(aes(color=n_entity), size=1.2) +
  scale_y_continuous(minor_breaks = seq(0 , 30, 2.5), breaks = seq(0, 30, 5)) +
  xlab("Situation") +
  ylab("Ratio RATP/RiRi")+
  labs(color="nombre d'entités")+
  theme(axis.text.x = element_text(angle = 90))


prow <- plot_grid(
  p1 + theme(legend.position = "none"), 
  p2+ theme(legend.position = "none"), 
  labels=c("A", "B"), 
  ncol = 2, nrow=1
)

legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0,0,0,12))
)

plegend <- plot_grid(prow, legend, rel_widths = c(3,.4))

title <- ggdraw() + draw_label("Temps de calcul de la lumière cumulé sur 180 jours de simulation", fontface='bold')

title

plot_grid(title, plegend,ncol=1, rel_heights = c(0.1, 1))



# courbe tot et run par steps
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/photomorpho_ramif_leg_sky46_1ent/')
df <- read.table("cputime_steps.csv", sep=',',stringsAsFactors = FALSE)
values <- as.numeric(df[[1]][2:length(df[[1]])])
steps <- 60:180
situation <- rep("RiRi l-egume", 121)

values <- append(values, as.numeric(df[[3]][2:length(df[[1]])]))
steps <- append(steps, 60:180)
situation <- append(situation, rep("RATP total", 121))

values <- append(values, as.numeric(df[[2]][2:length(df[[1]])]))
steps <- append(steps, 60:180)
situation <- append(situation, rep("RATP run", 121))

start <- 60
end <- 180
data_timesteps <- data.frame(x=sapply(steps,c), values=sapply(values,c), situation=sapply(situation, c))
# data_timesteps <- data_timesteps[data_timesteps$x > start & data_timesteps$x < end,]
p1 <- ggplot() + 
  geom_line(data=data_timesteps, aes(x=x, y=values, color=situation), size=1.2) +
  
  geom_ribbon(data=data_timesteps[data_timesteps$situation=="RATP total",], 
              aes(x = data_timesteps[data_timesteps$situation=="RATP total",1],
                  ymin=data_timesteps[data_timesteps$situation=="RATP run",2],
                  ymax=data_timesteps[data_timesteps$situation=="RATP total",2]),
                  alpha=0.2)+
  labs(x="Jour de l'année",y="CPU times en s")

  # geom_ribbon(data=data_timesteps[data_timesteps$situation=="RATP run",],
  #           aes(x = data_timesteps[data_timesteps$situation=="RATP run",1],
  #               ymin=data_timesteps[data_timesteps$situation=="RiRi l-egume",2],
  #               ymax=data_timesteps[data_timesteps$situation=="RATP run",2]),
  #               alpha=0.2)
  

start <- 70
end <- 85
data_timesteps_z1 <- data_timesteps[data_timesteps$x > start & data_timesteps$x < end,]
p2 <- ggplot() + 
  geom_line(data=data_timesteps_z1, aes(x=x, y=values, color=situation), size=1.2) +
  
  geom_ribbon(data=data_timesteps_z1[data_timesteps_z1$situation=="RATP total",], 
              aes(x = data_timesteps_z1[data_timesteps_z1$situation=="RATP total",1],
                  ymin=data_timesteps_z1[data_timesteps_z1$situation=="RATP run",2],
                  ymax=data_timesteps_z1[data_timesteps_z1$situation=="RATP total",2]),
              alpha=0.2)+
  labs(x="Jour de l'année",y="CPU times en s")

# geom_ribbon(data=data_timesteps[data_timesteps$situation=="RATP run",],
#             aes(x = data_timesteps[data_timesteps$situation=="RATP run",1],
#                 ymin=data_timesteps[data_timesteps$situation=="RiRi l-egume",2],
#                 ymax=data_timesteps[data_timesteps$situation=="RATP run",2]),
#             alpha=0.2)


start <- 120
end <- 135
data_timesteps_z2 <- data_timesteps[data_timesteps$x > start & data_timesteps$x < end,]
p3 <- ggplot() + 
  geom_line(data=data_timesteps_z2, aes(x=x, y=values, color=situation), size=1.2) +
  
  geom_ribbon(data=data_timesteps_z2[data_timesteps_z2$situation=="RATP total",], 
              aes(x = data_timesteps_z2[data_timesteps_z2$situation=="RATP total",1],
                  ymin=data_timesteps_z2[data_timesteps_z2$situation=="RATP run",2],
                  ymax=data_timesteps_z2[data_timesteps_z2$situation=="RATP total",2]),
              alpha=0.2)+
  labs(x="Jour de l'année",y="CPU times en s")

# geom_ribbon(data=data_timesteps[data_timesteps$situation=="RATP run",],
#             aes(x = data_timesteps[data_timesteps$situation=="RATP run",1],
#                 ymin=data_timesteps[data_timesteps$situation=="RiRi l-egume",2],
#                 ymax=data_timesteps[data_timesteps$situation=="RATP run",2]),
#             alpha=0.2)

# temps CPU en fonction du nombre de layers 
# tableau nb layer en fonction du step
layermax <- 0
steps <- list()
nlayers <- list()
for (n in 0:120)
{
  df <- read.table(paste(paste("diff_part", as.character(n), sep='_'), 'csv', sep='.'), sep=',', stringsAsFactors = FALSE)
  nlayers <- append(nlayers, ncol(df))
  steps <- append(steps, n)
}
df_nlayers_steps <- data.frame(doy=sapply(steps, c), nlayers=sapply(nlayers,c))

values <- list()
for (n in 0:120)
{
  values <- append(values, data_timesteps[data_timesteps$x==n+60 & data_timesteps$situation=="RATP run",2])
}
df_corr_cpu_layers <- data.frame(nlayers=sapply(nlayers,c), cputime=sapply(values,c))

df_line <- data.frame(x=c(0,40), y=c(0, max(df_corr_cpu_layers$cputime)))

p4 <- ggplot() +
  geom_line(data=df_line, aes(x=x, y=y), size=1.2, alpha=0.4, color="red") +
  geom_point(data=df_corr_cpu_layers, aes(x=nlayers, y=cputime), shape=3, size=2) +
  labs(x="Nombre de couches en z",y="CPU times en s", title="Analyse pour chaque jour de la simulation sur RATP run avec un ciel à 46 directions")
p4 + theme(plot.title = element_text(hjust=0.5))



legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0,0,0,12))
) + labs(colours="Fonction appelée")


pinit <- plot_grid(
  p1 + theme(legend.position = "none"), 
  p2 + theme(legend.position = "none"), 
  p3 + theme(legend.position = "none"), 
  legend,
  ncol=2, nrow=2, 
  labels=c('A', 'B', 'C'))

title <- ggdraw() + draw_label("Analyse du CPU time pour un ciel à 46 directions", fontface='bold')

title

plot_grid(title, pinit,ncol=1, rel_heights = c(0.1, 1))


## courbe PARa et PARt total l-egume etc...
# gauche para, droite part
file_list <- list(
  
  'photomorpho_ramif_leg_sky5_25t_1ent',
  'photomorpho_ramif_leg_sky46_1ent',
  'photomorpho_ramif_leg_sky100_1ent'
)
situation_list <- list("RATP: sky5 | 25 tubes", "RATP: sky46 | 25 tubes", "RATP: sky100 | 25 tubes")
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp')

values <- list()
situation <- list()
for (i in 1:3)
{
  df <- read.table(paste(file_list[[i]],'global_outputs.csv', sep='/'), sep=',',stringsAsFactors = FALSE)
  values <- append(values, as.numeric(format(round(as.numeric(df$V6[2]), 2), nsmall = 2)))
  situation <- append(situation, situation_list[[i]])
}
values <- append(values, as.numeric(format(round(as.numeric(df$V4[2]), 2), nsmall = 2)))
situation <- append(situation, "RiRi l-egume")

data_para <- data.frame(x=sapply(situation,c), values=sapply(values,c))
p1 <- ggplot(data=data_para, aes(x=reorder(x, values), y=values)) + 
  geom_col(fill="#0073C2FF", alpha=0.7) +
  geom_text(aes(label=values), vjust = 1.6, color="white")+
  xlab("Situation") +
  ylab("PARa en W/m²") +
  theme(axis.text.x = element_text(angle = 90))
p1 + scale_color_brewer(palette="Paired") + theme_minimal()

values <- list()
situation <- list()
for (i in 1:3)
{
  df <- read.table(paste(file_list[[i]],'global_outputs.csv', sep='/'), sep=',',stringsAsFactors = FALSE)
  values <- append(values, as.numeric(format(round(as.numeric(df$V7[2]), 2), nsmall = 2)))
  situation <- append(situation, situation_list[[i]])
}
values <- append(values, as.numeric(format(round(as.numeric(df$V5[2]), 2), nsmall = 2)))
situation <- append(situation, "RiRi l-egume")

data_part <- data.frame(x=sapply(situation,c), values=sapply(values,c))
p2 <- ggplot(data=data_part, aes(x=reorder(x, values), y=values)) + 
  geom_col(fill="#FC4E07", alpha=0.7) +
  geom_text(aes(label=values), vjust = 1.6, color="white")+
  xlab("Situation") +
  ylab("PARt en W/m²") +
  theme(axis.text.x = element_text(angle = 90))
p2 + scale_color_brewer(palette="Paired") + theme_minimal()

p_init <- plot_grid(p1, p2, ncol=2, nrow=1, labels=c("A","B"))

title <- ggdraw() + draw_label("Somme des rayonnements sur tous les jours de la simulation pour chaque modèle", fontface='bold')

title

plot_grid(title, p_init,ncol=1, rel_heights = c(0.1, 1))

## boxplot diff voxel par couche
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/photomorpho_ramif_leg_sky5_25t_1ent')

data_surf_layer <- plot_df_diffvoxel_perlayer('surf_vox_17112_l-egume', 120)

data_para_layer <- plot_df_diffvoxel_perlayer('diff_para_17112_l-egume', 120)
data_para_layer$rayonnement <- c(rep("PARa", length(data_para_layer[[1]])))

data_part_layer <- plot_df_diffvoxel_perlayer('diff_part', 120)
data_part_layer$rayonnement <- c(rep("PARt", length(data_part_layer[[1]])))

data_diff_vox_layers <- rbind(data_para_layer, data_part_layer)


p1 <- ggplot(data_diff_vox_layers, aes(group=rayonnement)) + 
      geom_point(aes(x,y, color=rayonnement),shape=18, size=2) +
      geom_line(aes(x, median, color=rayonnement)) +
      geom_ribbon(aes(x, ymin=min, ymax=max,fill=rayonnement, color=rayonnement), alpha=0.5) +
  labs(x="couche en z", y="Rayonnement (RATP - RiRi(l-egume)) en W/m²", title="Min, max et médiane des différence de rayonnement par couche sur tous les jours de simulation, RATP avec un ciel à 5 directions")

p2 <- ggplot(data_surf_layer) + 
  geom_point(aes(x,y),shape=18) +
  geom_line(aes(x, median)) +
  geom_ribbon(aes(x, ymin=min, ymax=max), alpha=0.4) +
  labs(x="couche en z", y="Surface foliaire en m²", title="Max et médiane des surfaces foliaires par couche sur tous les jours de simulation")

legend <- get_legend(p1  + theme(legend.box.margin = margin(0,0,0,5)))
prow <- plot_grid(
  p1 + theme(legend.position = "none"), 
  legend,
  p2,
  ncol=2, nrow=2,  aling='v', rel_widths = c(5, 1))
prow


## diff voxel en fonction de la surface ?

## boxplot diff voxel sur une couche par step
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/photomorpho_ramif_leg_sky5_25t_1ent')

layer <- 36

data_para_layer_steps <- plot_df_diffvoxel_persteps_layer('diff_para_17112_l-egume', 120, layer)
data_part_layer_steps <- plot_df_diffvoxel_persteps_layer('diff_part', 120, 36)
data_surf_layer_steps <- plot_df_diffvoxel_persteps_layer('surf_vox_17112_l-egume', 120, layer)

data_para_layer_steps$rayonnement <- c(rep("PARa", length(data_para_layer_steps[[1]])))
data_part_layer_steps$rayonnement <- c(rep("PARt", length(data_para_layer_steps[[1]])))

data_diff_vox_layers_steps <- rbind(data_para_layer_steps, data_part_layer_steps)

p1 <- ggplot(data_diff_vox_layers_steps, aes(group=rayonnement)) + 
  geom_point(aes(x,y, color=rayonnement),shape=18, size=2) +
  geom_line(aes(x, median, color=rayonnement)) +
  geom_ribbon(aes(x, ymin=min, ymax=max,fill=rayonnement, color=rayonnement), alpha=0.5) +
  labs(x="Jour de l'année", y="Rayonnement (RATP - RiRi(l-egume)) en W/m²", title="Min, max et médiane des différence de rayonnement sur la couche 36, évolution par jour, RATP avec un ciel à 5 directions")

p2 <- ggplot(data_surf_layer_steps) + 
  geom_point(aes(x,y),shape=18) +
  geom_line(aes(x, median)) +
  geom_ribbon(aes(x, ymin=min, ymax=max), alpha=0.4) +
  labs(x="Jour de l'année", y="Surface foliaire en m²", title="Max et médiane des surfaces foliaires sur la couche 36, évolution par jour")

legend <- get_legend(p1  + theme(legend.box.margin = margin(0,0,0,5)))
prow <- plot_grid(
  p1 + theme(legend.position = "none"), 
  legend,
  p2,
  ncol=2, nrow=2,  aling='v', rel_widths = c(5, 1))
prow

# diff en fonction des différences de rayonnement 
diff_max <- list()
diff_min <- list()
diff_mediane <- list()
rayonnement <- list()
surf_max <- list()
nsteps <- 120
for (i in 0:nsteps)
{
  if (as.numeric(i+60) %in%  data_para_layer_steps$x)
  {
    diff_max <- append(diff_max, data_para_layer_steps[data_para_layer_steps$x==(i+60),4][1])
    diff_min <- append(diff_min, data_para_layer_steps[data_para_layer_steps$x==(i+60),3][1])
    diff_mediane <- append(diff_mediane, data_para_layer_steps[data_para_layer_steps$x==(1+60),"median"])
    rayonnement <- append(rayonnement, "PARa")
    # diff_max <- append(diff_max, data_part_layer_steps[data_part_layer_steps$x==(i+60),4][1])
    # diff_min <- append(diff_min, data_part_layer_steps[data_part_layer_steps$x==(i+60),3][1])
    # diff_mediane <- append(diff_mediane, data_part_layer_steps[data_part_layer_steps$x==(i+60),"median"])
    # rayonnement <- append(rayonnement, "PARt")
    # surf_max <- append(surf_max, data_surf_layer_steps[data_surf_layer_steps$x==(i+60),4][1])
    surf_max <- append(surf_max, data_surf_layer_steps[data_surf_layer_steps$x==(i+60),4][1])
  }
}

df_dv_max <- data.frame(surface=sapply(surf_max, c), diff_rayon=sapply(diff_max, c), type=sapply(rayonnement, c))
df_line <- data.frame(x=c(0,max(df_dv_max$surface)), y=c(0, max(df_dv_max$diff_rayon)))

p1 <- ggplot() +
      geom_point(data=df_dv_max, aes(x=surface, y=diff_rayon, fill=type, color=type))+
      geom_line(data=df_line, aes(x=x, y=y), size=1.2, alpha=0.4, color="red")+
    labs(x="Surface foliaire en m²", y="Rayonnement (RiRi - RATP) en W/m²", title="Maximum des différences de rayonnement")
  
df_dv_min <- data.frame(surface=sapply(surf_max, c), diff_rayon=sapply(diff_min, c), type=sapply(rayonnement, c))
df_line <- data.frame(x=c(0,max(df_dv_min$surface)), y=c(0, min(df_dv_min$diff_rayon)))

p2 <- ggplot() +
  geom_point(data=df_dv_min, aes(x=surface, y=diff_rayon, fill=type, color=type))+
  geom_line(data=df_line, aes(x=x, y=y), size=1.2, alpha=0.4, color="red")+
  labs(x="Surface foliaire en m²", y="Rayonnement (RiRi - RATP) en W/m²", title="Minimum des différences de rayonnement")

plot_grid(p1, p2, ncol = 1)



# se place dans le dossier
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp')
legume_default_1 <- read.table('nophotomorpho_ramif_leg_sky5_1t/toto_17111_l-egume.csv', sep=';',stringsAsFactors = FALSE)
ratp_pass_sky5_1t_1 <- read.table('nophotomorpho_ramif_leg_sky5_1t/outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
ratp_pass_sky5_25t_1 <- read.table('nophotomorpho_ramif_leg_sky5_25t/outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
ratp_pass_sky46_1 <- read.table('nophotomorpho_ramif_leg_sky46_vtk/outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
ratp_pass_sky100_1 <- read.table('nophotomorpho_ramif_leg_sky100/outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)


diff <- read.table('nophotomorpho_ramif_leg_sky100/diff_para_17111_l-egume_0.csv', sep=',', stringsAsFactors = FALSE)


setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/photomorpho_leg_sky46_ratp_0/')
list_diff_voxels_para_0 <- list_dataframe_keyword("diff_para_17111_l-egume", 120)
list_diff_voxels_para_0 <- append(list_diff_voxels_para_0, list_dataframe_keyword("diff_para_17112_l-egume", 120))
list_diff_voxels_part_0 <- list_dataframe_keyword("diff_part", 120)

setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/photomorpho_leg_sky46_ratp_90/')
list_diff_voxels_para_90 <- list_dataframe_keyword("diff_para_17111_l-egume", 120)
list_diff_voxels_para_90 <- append(list_diff_voxels_para_0, list_dataframe_keyword("diff_para_17112_l-egume", 120))
list_diff_voxels_part_90 <- list_dataframe_keyword("diff_part", 120)

setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/photomorpho_leg_sky46_ratp_-90/')
list_diff_voxels_para_m90 <- list_dataframe_keyword("diff_para_17111_l-egume", 120)
list_diff_voxels_para_m90 <- append(list_diff_voxels_para_0, list_dataframe_keyword("diff_para_17112_l-egume", 120))
list_diff_voxels_part_m90 <- list_dataframe_keyword("diff_part", 120)

setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/photomorpho_leg_sky46_ratp_180/')
list_diff_voxels_para_180 <- list_dataframe_keyword("diff_para_17111_l-egume", 120)
list_diff_voxels_para_180 <- append(list_diff_voxels_para_0, list_dataframe_keyword("diff_para_17112_l-egume", 120))
list_diff_voxels_part_180 <- list_dataframe_keyword("diff_part", 120)

df_diffvox_para_allsteps <- plot_dataframe_hist_alltimesteps(list(list_diff_voxels_para_0, 
                                                            list_diff_voxels_para_90, 
                                                            list_diff_voxels_para_m90, 
                                                            list_diff_voxels_para_180), list("0", "+pi/2", "-pi/2", "pi"))
p <- ggplot(df_diffvox_para_allsteps, aes(x=x, y=values)) + geom_boxplot() + 
  labs(title="PARa (RiRi - RATP) sur tous les pas de temps", x="situation", y="[0, 1]")
p

df_diffvox_part_allsteps <- plot_dataframe_hist_alltimesteps(list(list_diff_voxels_part_0, 
                                                                  list_diff_voxels_part_90, 
                                                                  list_diff_voxels_part_m90, 
                                                                  list_diff_voxels_part_180), list("0", "+pi/2", "-pi/2", "pi"))
p <- ggplot(df_diffvox_part_allsteps, aes(x=x, y=values)) + geom_boxplot() + 
  labs(title="PARt (RiRi - RATP) sur tous les pas de temps", x="situation", y="[0, 1]")
p


setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp/photomorpho_leg_sky46_ratp_90/')

i <- 3
values <- list()
filename <- paste(paste("diff_para_17111_l-egume", as.character(i), sep="_"), "csv", sep=".")
df <- read.table(filename, sep=',', stringsAsFactors = FALSE)
values <- append(values, sapply(df[2:length(df),], function(x) as.numeric(as.character(x))))


df <- plot_dataframe_diffvoxel_perstep("diff_part", 120, list(max, min), list("max", "min"))
p <- ggplot(data=df, aes(x=x, y=values, group=names)) + 
  geom_line(aes(color=names)) +
  ggtitle('Epsi sommÃ© sur le couvert') +
  xlab("DOY") +
  ylab("epsi")
p


## Epsi ##
# df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp_passive), list("default", "ratp"))
# df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp_passive), list("default", "ratp"))
# df <- plot_dataframe_canopy("epsi", list(legume_ratp_passive), list("ratp"))
df <- plot_dataframe_canopy("epsi", list(legume_default_1, 
                                         ratp_pass_sky5_1t_1, 
                                         ratp_pass_sky5_25t_1, 
                                         ratp_pass_sky46_1,
                                         ratp_pass_sky100_1), list("default", 
                                                                      "ratp passive sky_5 1 tube", 
                                                                      "ratp passive sky_5 25 tube",
                                                                      "ratp passive sky_46 25 tube",
                                                                      "ratp passive sky_100 25 tube"), 175, 178)

p <- ggplot(data=df, aes(x=x, y=values, group=light)) + 
      geom_line(aes(color=light)) +
      ggtitle('Epsi sommÃ© sur le couvert') +
      xlab("DOY") +
      ylab("epsi")
p

df <- plot_dataframe_plant("epsi", list(legume_default_1, 
                                         ratp_pass_sky5_1t_1), list("default", 
                                                                   "ratp passive sky_5 1 tube"), list(10, 12), 100, 178)

p <- ggplot(data=df, aes(x=x, y=values, group=light)) + 
  geom_line(aes(color=light)) +
  ggtitle('Epsi') +
  xlab("DOY") +
  ylab("epsi")
p

## SurfPlante ##
# df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp_passive), list("default", "ratp"))
df <- plot_dataframe_canopy("SurfPlante", list(legume_default, legume_ratp_s5, legume_ratp_s46), list("default", "ratp active ciel sky_5", "ratp active ciel turtle46"))

p <- ggplot(data=df, aes(x=x, y=values, group=light)) + 
  geom_line(aes(color=light)) +
  ggtitle('SurfPlante sommÃ© sur le couvert') +
  xlab("DOY") +
  ylab("Surface foliaire en mÂ²")
p

## Hplante ##
df <- plot_dataframe_canopy("Hplante", list(legume_default, legume_ratp_s5, legume_ratp_s46), list("default", "ratp active ciel sky_5", "ratp active ciel turtle46"))

p <- ggplot(data=df, aes(x=x, y=values, group=light)) + 
  geom_line(aes(color=light)) +
  ggtitle('Hplante sommÃ© sur le couvert') +
  xlab("DOY") +
  ylab("Hauteur en m")
p


## COUVERT ##
setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/legume_ratp')
# 4 graphes par situation avec toutes les situations
# ratp passive
photo_ram_legume_default_1 <- read.table('photomorpho_ramif_leg_sky5_1t/toto_17111_l-egume.csv', sep=';',stringsAsFactors = FALSE)
pass_ratp_sky5_1t_1 <- read.table('photomorpho_ramif_leg_sky5_1t/outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
pass_ratp_sky5_25t_1 <- read.table('photomorpho_ramif_leg_sky5_25t/outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
pass_ratp_sky46_1 <- read.table('photomorpho_ramif_leg_sky46_vtk/outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)
pass_ratp_sky100_1 <- read.table('photomorpho_ramif_leg_sky100/outputs_ratp_passive_17111_l-egume.csv', sep=',', stringsAsFactors = FALSE)

diff_pass_ratp_sky5_1t_1 <- diff_dataframe_toto_legume(photo_ram_legume_default_1, pass_ratp_sky5_1t_1, passive=TRUE)
diff_pass_ratp_sky5_25t_1 <- diff_dataframe_toto_legume(photo_ram_legume_default_1, pass_ratp_sky5_25t_1, passive=TRUE)
diff_pass_ratp_sky46_1 <- diff_dataframe_toto_legume(photo_ram_legume_default_1, pass_ratp_sky46_1, passive=TRUE)
diff_pass_ratp_sky100_1 <- diff_dataframe_toto_legume(photo_ram_legume_default_1, pass_ratp_sky100_1, passive=TRUE)

photo_ram_legume_default_2 <- read.table('photomorpho_ramif_leg_sky5_1t/toto_17112_l-egume.csv', sep=';',stringsAsFactors = FALSE)
pass_ratp_sky5_1t_2 <- read.table('photomorpho_ramif_leg_sky5_1t/outputs_ratp_passive_17112_l-egume.csv', sep=',', stringsAsFactors = FALSE)
pass_ratp_sky5_25t_2 <- read.table('photomorpho_ramif_leg_sky5_25t/outputs_ratp_passive_17112_l-egume.csv', sep=',', stringsAsFactors = FALSE)
pass_ratp_sky46_2 <- read.table('photomorpho_ramif_leg_sky46_vtk/outputs_ratp_passive_17112_l-egume.csv', sep=',', stringsAsFactors = FALSE)
pass_ratp_sky100_2 <- read.table('photomorpho_ramif_leg_sky100/outputs_ratp_passive_17112_l-egume.csv', sep=',', stringsAsFactors = FALSE)

diff_pass_ratp_sky5_1t_2 <- diff_dataframe_toto_legume(photo_ram_legume_default_2, pass_ratp_sky5_1t_2, passive=TRUE)
diff_pass_ratp_sky5_25t_2 <- diff_dataframe_toto_legume(photo_ram_legume_default_2, pass_ratp_sky5_25t_2, passive=TRUE)
diff_pass_ratp_sky46_2 <- diff_dataframe_toto_legume(photo_ram_legume_default_2, pass_ratp_sky46_2, passive=TRUE)
diff_pass_ratp_sky100_2 <- diff_dataframe_toto_legume(photo_ram_legume_default_2, pass_ratp_sky100_2, passive=TRUE)

# ratp active photomorpho ramif
photo_ram_ratp_sky5_1t_1 <- read.table('photomorpho_ramif_ratp_sky5_1t/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_ram_ratp_sky5_25t_1 <- read.table('photomorpho_ramif_ratp_sky5_25t/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_ram_ratp_sky46_1 <- read.table('photomorpho_ramif_ratp_sky46_vtk/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_ram_ratp_sky100_1 <- read.table('photomorpho_ramif_ratp_sky100/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)

diff_photo_ram_ratp_sky5_1t_1 <- diff_dataframe_toto_legume(photo_ram_legume_default_1, photo_ram_ratp_sky5_1t_1)
diff_photo_ram_ratp_sky5_25t_1 <- diff_dataframe_toto_legume(photo_ram_legume_default_1, photo_ram_ratp_sky5_25t_1)
diff_photo_ram_ratp_sky46_1 <- diff_dataframe_toto_legume(photo_ram_legume_default_1, photo_ram_ratp_sky46_1)
diff_photo_ram_ratp_sky100_1 <- diff_dataframe_toto_legume(photo_ram_legume_default_1, photo_ram_ratp_sky100_1)

photo_ram_ratp_sky5_1t_2 <- read.table('photomorpho_ramif_ratp_sky5_1t/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_ram_ratp_sky5_25t_2 <- read.table('photomorpho_ramif_ratp_sky5_25t/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_ram_ratp_sky46_2 <- read.table('photomorpho_ramif_ratp_sky46_vtk/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_ram_ratp_sky100_2 <- read.table('photomorpho_ramif_ratp_sky100/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)

diff_photo_ram_ratp_sky5_1t_2 <- diff_dataframe_toto_legume(photo_ram_legume_default_2, photo_ram_ratp_sky5_1t_2)
diff_photo_ram_ratp_sky5_25t_2 <- diff_dataframe_toto_legume(photo_ram_legume_default_2, photo_ram_ratp_sky5_25t_2)
diff_photo_ram_ratp_sky46_2 <- diff_dataframe_toto_legume(photo_ram_legume_default_2, photo_ram_ratp_sky46_2)
diff_photo_ram_ratp_sky100_2 <- diff_dataframe_toto_legume(photo_ram_legume_default_2, photo_ram_ratp_sky100_2)

# ratp active photomorpho noramif
photo_noram_legume_default_1 <- read.table('photomorpho_noramif_legume_vtk/toto_17111_l-egume.csv', sep=';',stringsAsFactors = FALSE)
photo_noram_ratp_sky5_1t_1 <- read.table('photomorpho_noramif_ratp_sky5_1t/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_noram_ratp_sky5_25t_1 <- read.table('photomorpho_noramif_ratp_sky5_25t/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_noram_ratp_sky46_1 <- read.table('photomorpho_noramif_ratp_sky46_vtk/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_noram_ratp_sky100_1 <- read.table('photomorpho_noramif_ratp_sky100/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)

diff_photo_noram_ratp_sky5_1t_1 <- diff_dataframe_toto_legume(photo_noram_legume_default_1, photo_noram_ratp_sky5_1t_1)
diff_photo_noram_ratp_sky5_25t_1 <- diff_dataframe_toto_legume(photo_noram_legume_default_1, photo_noram_ratp_sky5_25t_1)
diff_photo_noram_ratp_sky46_1 <- diff_dataframe_toto_legume(photo_noram_legume_default_1, photo_noram_ratp_sky46_1)
diff_photo_noram_ratp_sky100_1 <- diff_dataframe_toto_legume(photo_noram_legume_default_1, photo_noram_ratp_sky100_1)

photo_noram_legume_default_2 <- read.table('photomorpho_noramif_legume_vtk/toto_17112_l-egume.csv', sep=';',stringsAsFactors = FALSE)
photo_noram_ratp_sky5_1t_2 <- read.table('photomorpho_noramif_ratp_sky5_1t/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_noram_ratp_sky5_25t_2 <- read.table('photomorpho_noramif_ratp_sky5_25t/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_noram_ratp_sky46_2 <- read.table('photomorpho_noramif_ratp_sky46_vtk/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
photo_noram_ratp_sky100_2 <- read.table('photomorpho_noramif_ratp_sky100/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)

diff_photo_noram_ratp_sky5_1t_2 <- diff_dataframe_toto_legume(photo_noram_legume_default_2, photo_noram_ratp_sky5_1t_2)
diff_photo_noram_ratp_sky5_25t_2 <- diff_dataframe_toto_legume(photo_noram_legume_default_2, photo_noram_ratp_sky5_25t_2)
diff_photo_noram_ratp_sky46_2 <- diff_dataframe_toto_legume(photo_noram_legume_default_2, photo_noram_ratp_sky46_2)
diff_photo_noram_ratp_sky100_2 <- diff_dataframe_toto_legume(photo_noram_legume_default_2, photo_noram_ratp_sky100_2)

# ratp active nophotomorpho
nophoto_ram_legume_default_1 <- read.table('nophotomorpho_ramif_leg_sky5_1t/toto_17111_l-egume.csv', sep=';',stringsAsFactors = FALSE)
nophoto_ram_ratp_sky5_1t_1 <- read.table('nophotomorpho_ramif_ratp_sky5_1t/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
nophoto_ram_ratp_sky5_25t_1 <- read.table('nophotomorpho_ramif_ratp_sky5_25t/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
nophoto_ram_ratp_sky46_1 <- read.table('nophotomorpho_ramif_ratp_sky46_vtk/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)
nophoto_ram_ratp_sky100_1 <- read.table('nophotomorpho_ramif_ratp_sky100/toto_17111_l-egume.csv', sep=';', stringsAsFactors = FALSE)

diff_nophoto_ram_ratp_sky5_1t_1 <- diff_dataframe_toto_legume(nophoto_ram_legume_default_1, nophoto_ram_ratp_sky5_1t_1)
diff_nophoto_ram_ratp_sky5_25t_1 <- diff_dataframe_toto_legume(nophoto_ram_legume_default_1, nophoto_ram_ratp_sky5_25t_1)
diff_nophoto_ram_ratp_sky46_1 <- diff_dataframe_toto_legume(nophoto_ram_legume_default_1, nophoto_ram_ratp_sky46_1)
diff_nophoto_ram_ratp_sky100_1 <- diff_dataframe_toto_legume(nophoto_ram_legume_default_1, nophoto_ram_ratp_sky100_1)

nophoto_ram_legume_default_2 <- read.table('nophotomorpho_ramif_leg_sky5_1t/toto_17112_l-egume.csv', sep=';',stringsAsFactors = FALSE)
nophoto_ram_ratp_sky5_1t_2 <- read.table('nophotomorpho_ramif_ratp_sky5_1t/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
nophoto_ram_ratp_sky5_25t_2 <- read.table('nophotomorpho_ramif_ratp_sky5_25t/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
nophoto_ram_ratp_sky46_2 <- read.table('nophotomorpho_ramif_ratp_sky46_vtk/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)
nophoto_ram_ratp_sky100_2 <- read.table('nophotomorpho_ramif_ratp_sky100/toto_17112_l-egume.csv', sep=';', stringsAsFactors = FALSE)

diff_nophoto_ram_ratp_sky5_1t_2 <- diff_dataframe_toto_legume(nophoto_ram_legume_default_2, nophoto_ram_ratp_sky5_1t_2)
diff_nophoto_ram_ratp_sky5_25t_2 <- diff_dataframe_toto_legume(nophoto_ram_legume_default_2, nophoto_ram_ratp_sky5_25t_2)
diff_nophoto_ram_ratp_sky46_2 <- diff_dataframe_toto_legume(nophoto_ram_legume_default_2, nophoto_ram_ratp_sky46_2)
diff_nophoto_ram_ratp_sky100_2 <- diff_dataframe_toto_legume(nophoto_ram_legume_default_2, nophoto_ram_ratp_sky100_2)

# Listes des dataframes (données bruts)
c_listes_entite_1 <- list(
  list(photo_ram_legume_default_1, 
       pass_ratp_sky5_1t_1, 
       pass_ratp_sky5_25t_1, 
       pass_ratp_sky46_1,
       pass_ratp_sky100_1),
  list(photo_noram_legume_default_1, 
       photo_noram_ratp_sky5_1t_1, 
       photo_noram_ratp_sky5_25t_1, 
       photo_noram_ratp_sky46_1,
       photo_noram_ratp_sky100_1),
  list(photo_ram_legume_default_1, 
       photo_ram_ratp_sky5_1t_1, 
       photo_ram_ratp_sky5_25t_1, 
       photo_ram_ratp_sky46_1,
       photo_ram_ratp_sky100_1),
  list(nophoto_ram_legume_default_1, 
       nophoto_ram_ratp_sky5_1t_1, 
       nophoto_ram_ratp_sky5_25t_1, 
       nophoto_ram_ratp_sky46_1,
       nophoto_ram_ratp_sky100_1)
)

c_listes_entite_2 <- list(
  list(photo_ram_legume_default_2, 
       pass_ratp_sky5_1t_2, 
       pass_ratp_sky5_25t_2, 
       pass_ratp_sky46_2,
       pass_ratp_sky100_2),
  list(photo_noram_legume_default_2, 
       photo_noram_ratp_sky5_1t_2, 
       photo_noram_ratp_sky5_25t_2, 
       photo_noram_ratp_sky46_2,
       photo_noram_ratp_sky100_2),
  list(photo_ram_legume_default_2, 
       photo_ram_ratp_sky5_1t_2, 
       photo_ram_ratp_sky5_25t_2, 
       photo_ram_ratp_sky46_2,
       photo_ram_ratp_sky100_2),
  list(nophoto_ram_legume_default_2, 
       nophoto_ram_ratp_sky5_1t_2, 
       nophoto_ram_ratp_sky5_25t_2, 
       nophoto_ram_ratp_sky46_2,
       nophoto_ram_ratp_sky100_2)
)

liste_legend_name <- list("default", 
                          "ratp sky5 1 tube", 
                          "ratp sky5 25 tube",
                          "ratp sky46 25 tube",
                          "ratp sky100 25 tube")

#### epsi
varname <- "epsi"
ytitle <- "Epsi"
globaltitle <- "Somme du epsilon sur tout le couvert, ZOOM"
start <- 160
end <- 180
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=TRUE)
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

#### Hauteur
varname <- "Hplante"
ytitle <- "Hauteur en m"

start <- 60
end <- 180
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=FALSE)
p_grid

start <- 160
globaltitle <- "Somme des hauteurs des plantes sur tout le couvert, ZOOM"
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=FALSE)
p_grid

#### MStot (demande en azote des parties aériennes)
varname <- "MStot"
ytitle <- "Demande en azote en T.ha-1"
globaltitle <- "Somme des demandes en azote des parties aériennes sur tout le couvert"
start <- 60
end <- 180
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=FALSE)
p_grid

globaltitle <- "Somme des demandes en azote des parties aériennes sur tout le couvert, ZOOM"
start <- 160
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=FALSE)
p_grid

#### PARi
varname <- "PARiPlante"
ytitle <- "PARi en W.m-2"
globaltitle <- "Somme du PAR intercepté  sur tout le couvert"
start <- 60
end <- 180
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=TRUE)
p_grid

globaltitle <- "Somme du PAR intercepté  sur tout le couvert, ZOOM"
start <- 160
p_grid <- plot_variable_canopy(varname, start, end, c_listes_entite_1, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=TRUE)
p_grid

#### Nouvelles listes : différences (l-egume - ratp)
# Listes des dataframes (données bruts)
c_listes_entite_1 <- list(
  list(diff_pass_ratp_sky5_1t_1, 
       diff_pass_ratp_sky5_25t_1, 
       diff_pass_ratp_sky46_1,
       diff_pass_ratp_sky100_1),
  list(diff_photo_noram_ratp_sky5_1t_1, 
       diff_photo_noram_ratp_sky5_25t_1, 
       diff_photo_noram_ratp_sky46_1,
       diff_photo_noram_ratp_sky100_1),
  list(diff_photo_ram_ratp_sky5_1t_1, 
       diff_photo_ram_ratp_sky5_25t_1, 
       diff_photo_ram_ratp_sky46_1,
       diff_photo_ram_ratp_sky100_1),
  list(diff_nophoto_ram_ratp_sky5_1t_1, 
       diff_nophoto_ram_ratp_sky5_25t_1, 
       diff_nophoto_ram_ratp_sky46_1,
       diff_nophoto_ram_ratp_sky100_1)
)

c_listes_entite_2 <- list(
  list(diff_pass_ratp_sky5_1t_2, 
       diff_pass_ratp_sky5_25t_2, 
       diff_pass_ratp_sky46_2,
       diff_pass_ratp_sky100_2),
  list(diff_photo_noram_ratp_sky5_1t_2, 
       diff_photo_noram_ratp_sky5_25t_2, 
       diff_photo_noram_ratp_sky46_2,
       diff_photo_noram_ratp_sky100_2),
  list(diff_photo_ram_ratp_sky5_1t_2, 
       diff_photo_ram_ratp_sky5_25t_2, 
       diff_photo_ram_ratp_sky46_2,
       diff_photo_ram_ratp_sky100_2),
  list(diff_nophoto_ram_ratp_sky5_1t_2, 
       diff_nophoto_ram_ratp_sky5_25t_2, 
       diff_nophoto_ram_ratp_sky46_2,
       diff_nophoto_ram_ratp_sky100_2)
)

liste_legend_name <- list("ratp sky5 1 tube", 
                          "ratp sky5 25 tube",
                          "ratp sky46 25 tube",
                          "ratp sky100 25 tube")
#### PAR ENTITE
# difference legume - ratp par espèce
# LAI par espèce (pour voir si il y a pas une espèce qui déconne plus que l'autre)
varname <- "SurfPlante"
ytitle <- "Différence (l-egume - ratp) : LAI"
globaltitle <- "Différence (l-egume - ratp) sur le LAI sur l'entité 2"
start <- 60
end <- 180
surfsol <- 0.08
p_grid <- plot_variable_entity(varname, start, end, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=FALSE, lai_surfsol = surfsol)
p_grid

# epsi
varname <- "epsi"
ytitle <- "Différence (l-egume - ratp) : Epsi"
globaltitle <- "Différence (l-egume - ratp) sur le Epsilon sur l'entité 2"
start <- 60
end <- 180

p_grid <- plot_variable_entity(varname, start, end, c_listes_entite_2, liste_legend_name, ytitle, globaltitle, passive=TRUE)
p_grid


#### PAR PLANTE
# pour chaque variable
#   
# stats : 
#   pour chaque plante :
#       pour chaque situation
#         legume - ratp

# liste de plantes, sur un ciel ratp particulier
# étude sur RATP ciel à 46 directions
c_listes_entite_1 <- list(
  list(diff_pass_ratp_sky46_1),
  list(diff_photo_noram_ratp_sky46_1),
  list(diff_photo_ram_ratp_sky46_1),
  list(diff_nophoto_ram_ratp_sky46_1)
)

c_listes_entite_2 <- list(
  list(diff_pass_ratp_sky46_2),
  list(diff_photo_noram_ratp_sky46_2),
  list(diff_photo_ram_ratp_sky46_2),
  list(diff_nophoto_ram_ratp_sky46_2)
)

# PARi
varname <- "PARiPlante"
ytitle <- "PAR en W.m-2"
globaltitle <- "Différence (l-egume - ratp) sur le PAR intercepté, entité 1"
start <- 60
end <- 180
list_legend <- list("ratp sky46")

## entité 1
# 02
list_plantes <- list(3,4, 5, 6, 7, 8)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_1, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# 06-09
list_plantes <- list(9, 10, 11, 12, 13, 14, 15)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_1, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# 014-018-017
list_plantes <- list(16,17, 18, 19, 20, 21)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_1, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# 019-023-024
list_plantes <- list(22, 23, 24, 25, 26, 27, 28)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_1, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# 030-027-28
list_plantes <- list(29,30, 31, 32, 33, 34)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_1, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# max
list_plantes <- list(5, 26,20,27,17)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_1, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# 023-017-019-02-027-024-014-

## entite 2
# 02
list_plantes <- list(3,4, 5, 6, 7, 8)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_2, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# 08
list_plantes <- list(9, 10, 11, 12, 13, 14, 15)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_2, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

#  014-015
list_plantes <- list(16,17, 18, 19, 20, 21)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_2, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# 020-024
list_plantes <- list(22, 23, 24, 25, 26, 27, 28)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_2, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# 026-027
list_plantes <- list(29,30, 31, 32, 33, 34)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_2, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# max
globaltitle <- "Différence (l-egume - ratp) sur le PAR intercepté, entité 2"
list_plantes <- list(5, 17, 27, 30)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_2, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

# epsi
varname <- "epsi"
ytitle <- "Epsi"
globaltitle <- "Différence (l-egume - ratp) sur Epsilon, entité 1"
start <- 60
end <- 180
list_legend <- list("ratp sky46")

list_plantes <- list(5, 26,20,27,17)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_1, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid

globaltitle <- "Différence (l-egume - ratp) sur Epsilon, entité 2"
list_plantes <- list(5, 17, 27, 30)
p_grid <- plot_variable_plant(varname, start, end, c_listes_entite_1, list_legend, list_plantes, ytitle, globaltitle, passive=TRUE)
p_grid


# voir max hplante
varname <- "Hplante"
df <- max_per_plant(var_name, c_listes_entite_2[[1]][[1]])
print("default, photo ramif")
print(df[order(df$max),])

print("default, photo, noram")
df <- max_per_plant(var_name, c_listes_entite_2[[2]][[1]])
print(df[order(df$max),])

print("default, nophoto")
df <- max_per_plant(var_name, c_listes_entite_2[[4]][[1]])
print(df[order(df$max),])

print("ratp sky46, photo ram")
df <- max_per_plant(var_name, c_listes_entite_2[[3]][[4]])
print(df[order(df$max),])


