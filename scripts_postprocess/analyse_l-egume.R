setwd('C:/Users/mwoussen/cdd/codes/vegecouplelight/scripts_postprocess')
source("my_functions.R")

library(ggplot2)



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
situation_list <- list("RATP: sky5 | 1tube", "RATP: sky5 | 25 tubes", "RATP: sky46 | 25 tubes", "RATP: sky100 | 25 tubes")
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
situation <- append(situation, "l-egume")
nb_entity <- append(nb_entity, "2 entités") 

df <- read.table(file_list[[5]], sep=',',stringsAsFactors = FALSE)
values <- append(values, as.numeric(df$V1[2]))
situation <- append(situation, "l-egume")
nb_entity <- append(nb_entity, "1 entité") 

data_cputime <- data.frame(x=sapply(situation,c), values=sapply(values,c), n_entity=sapply(nb_entity, c))

p <- ggplot(data=data_cputime, aes(x=reorder(x, values), y=values, group=n_entity)) + 
  geom_line(aes(color=n_entity)) +
  ggtitle('Epsi sommÃ© sur le couvert') +
  xlab("Situation") +
  ylab("CPU time en s")
p
























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



      
