### TRAITEMENT SORTIES L-EGUME ###


########################
#                      #
#     A NETTOYER       #
#                      #
########################



# courbes somme sur le couvert
# inputs :
#   - liste de dataframe
#   - nom variable Ã  comparer
#   - prÃ©ciser si temps rÃ©duit
plot_dataframe_canopy <- function(variable, dataframe_list, name_list,t_start=-1, t_end=-1)
{
  # nouvelle liste colonne : steps, value, name
  values <- list()
  steps <- list()
  names <- list()
  for (i in 1:length(dataframe_list))
  {
    # résultat par plante
    if (ncol(dataframe_list[[i]]) > 5)
    {
      values <- append(values, rowSums(sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,3:ncol(dataframe_list[[i]])], function(x) as.numeric(as.character(x)))))
    }
    # tableau avec juste le résultat sur le couvert
    else
    {
      values <- append(values, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,3], function(x) as.numeric(as.character(x))))
    }
    steps <- append(steps, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,2], function(x) as.numeric(as.character(x))))
    names <- append(names, rep(name_list[[i]], length(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,2])))
  }


  df <- data.frame(x=sapply(steps,c), values=sapply(values,c), light=sapply(names,c))
  if (t_start > -1)
  {
    if (t_end > -1)
    {
      df <- df[df$x >= t_start & df$x <= t_end,]
    }
    df <- df[df$x >= t_start,]
  }

  df
}

# affiche par plante : list de num de plantes
plot_dataframe_plant <- function(variable, dataframe_list, name_list, plant_list, t_start=-1, t_end=-1)
{
  # nouvelle liste colonne : steps, value, name
  values <- list()
  steps <- list()
  names <- list()
  for (i in 1:length(dataframe_list))
  {
    for (j in 1:length(plant_list))
    {
      values <- append(values, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,plant_list[[j]]], function(x) as.numeric(as.character(x))))
      steps <- append(steps, sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,2], function(x) as.numeric(as.character(x))))
      names <- append(names, rep(name_list[[i]], #paste(name_list[[i]], as.character(plant_list[[j]] -3), sep=" plante "),
                                 length(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,2])))
    }
  }
  df <- data.frame(x=sapply(steps,c), values=sapply(values,c), light=sapply(names,c))
  if (t_start > -1)
  {
    if (t_end > -1)
    {
      df <- df[df$x >= t_start & df$x <= t_end,]
    }
    df <- df[df$x >= t_start,]
  }

  df
}

# analyse diff des voxels sur tous les pas de temps
plot_dataframe_hist_alltimesteps <- function(list_list_list, situation_names){
  values <- list()
  steps <- list()
  for (k in 1:length(situation_names))
  {
    dataframe_list_para <- list_list_list[[k]]
    for (i in 1:length(dataframe_list_para))
    {
      values <- append(values, sapply(dataframe_list_para[[i]][2:length(dataframe_list_para[[i]])], function(x) as.numeric(as.character(x))))
      steps <- append(steps, rep(situation_names[[k]], length(dataframe_list_para[[i]])-1))
    }
  }

  df <- data.frame(x=sapply(steps,c), values=sapply(values,c))
  df
}

# crée une liste de listes à partir d'un mot clé sur des steps
list_dataframe_keyword <- function(keyword, nsteps){
  df_list <- list()
  for (i in 0:nsteps)
  {
    filename <- paste(paste(keyword, as.character(i), sep="_"), "csv", sep=".")
    df_list <- append(df_list, read.table(filename, sep=',', stringsAsFactors = FALSE))
  }
  df_list
}

# courbe ou hist ?
plot_dataframe_diffvoxel_perstep <- function(keyword, nsteps, functions, func_names){
  values <- list()
  steps <- list()
  names <- list()
  for (k in 1:length(functions))
  {
    for (i in 0:nsteps)
    {
      filename <- paste(paste(keyword, as.character(i), sep="_"), "csv", sep=".")
      df <- read.table(filename, sep=',', stringsAsFactors = FALSE)
      modify_input <- sapply(df[[1]][2:nrow(df)], function(x) as.numeric(as.character(x)))
      values <- append(values, functions[[k]](modify_input))
      steps <- append(steps, as.numeric(i))
      names <- append(names, func_names[[k]])

    }
  }

  df_final <- data.frame(x=sapply(steps,c), values=sapply(values,c), names=sapply(names, c))
  df_final
}

# on est déjà dans le dossier
plot_df_diffvoxel_perlayer <- function(keyword, nsteps){
  values <- list()
  layer <- list()
  layermax <- 0

  # trouver layermax
  for (n in 0:nsteps)
  {
    df <- read.table(paste(paste(keyword, as.character(n), sep='_'), 'csv', sep='.'), sep=',', stringsAsFactors = FALSE)
    if (ncol(df) > layermax) {layermax <- ncol(df)}
  }

  for (n in 0:nsteps)
  {
    df <- read.table(paste(paste(keyword, as.character(n), sep='_'), 'csv', sep='.'), sep=',', stringsAsFactors = FALSE)

    for (i in 1:ncol(df))
    {
      values <- append(values, as.numeric(df[[i]][2:length(df[[i]])]))
      layer <- append(layer, rep(as.character(layermax - (ncol(df) -i)  ), 100))
    }
  }

  df <- data.frame(x=sapply(layer, c), y=sapply(values,c))

  # reprend juste min-max et quelques valeurs entre
  values <- list()
  min <- list()
  max <- list()
  median <- list()
  layer <- list()
  for (i in 1:layermax)
  {
    values <- append(values, max(df[df$x==as.character(i), 2]))
    values <- append(values, min(df[df$x==as.character(i), 2]))
    values <- append(values, mean(df[df$x==as.character(i), 2]))
    values <- append(values, median(df[df$x==as.character(i), 2]))
    values <- append(values, mean(df[df$x==as.character(i), 2]) + var(df[df$x==as.character(i), 2]))
    values <- append(values, mean(df[df$x==as.character(i), 2]) - var(df[df$x==as.character(i), 2]))

    min <- append(min, rep(min(df[df$x==as.character(i), 2]), 6))
    max <- append(max, rep(max(df[df$x==as.character(i), 2]), 6))
    median <- append(median, rep(median(df[df$x==as.character(i), 2]), 6))


    layer <- append(layer, rep(i, 6))
  }

  df_final <- data.frame(x=sapply(layer, c), y=sapply(values,c), min=sapply(min, c), max=sapply(max, c), median=sapply(median, c))
  df_final
}

# on est déjà dans le dossier
plot_df_diffvoxel_persteps_layer <- function(keyword, nsteps, layer){

  # trouver layermax
  layermax <- 0
  for (n in 0:nsteps)
  {
    df <- read.table(paste(paste(keyword, as.character(n), sep='_'), 'csv', sep='.'), sep=',', stringsAsFactors = FALSE)
    if (ncol(df) > layermax) {layermax <- ncol(df)}
  }

  values <- list()
  steps <- list()
  nstart <- nsteps+1
  for (n in 0:nsteps)
  {
    df <- read.table(paste(paste(keyword, as.character(n), sep='_'), 'csv', sep='.'), sep=',', stringsAsFactors = FALSE)

    layer_rel <- layer + ncol(df) - layermax

    if (layer_rel > 0)
    {
      values <- append(values, as.numeric(df[[layer_rel]][2:length(df[[layer_rel]])]))
      steps <- append(steps, rep(as.character(n), 100))
      nstart <- pmin(n, nstart)
    }
  }

  df <- data.frame(x=sapply(steps, c), y=sapply(values,c))
  df

  # reprend juste min-max et quelques valeurs entre
  values <- list()
  min <- list()
  max <- list()
  median <- list()
  steps <- list()
  for (i in nstart:nsteps)
  {
    values <- append(values, max(df[df$x==as.character(i), 2]))
    values <- append(values, min(df[df$x==as.character(i), 2]))
    values <- append(values, mean(df[df$x==as.character(i), 2]))
    values <- append(values, median(df[df$x==as.character(i), 2]))
    values <- append(values, mean(df[df$x==as.character(i), 2]) + var(df[df$x==as.character(i), 2]))
    values <- append(values, mean(df[df$x==as.character(i), 2]) - var(df[df$x==as.character(i), 2]))

    min <- append(min, rep(min(df[df$x==as.character(i), 2]), 6))
    max <- append(max, rep(max(df[df$x==as.character(i), 2]), 6))
    median <- append(median, rep(median(df[df$x==as.character(i), 2]), 6))

    steps <- append(steps, rep(i+60, 6))
  }

  df_final <- data.frame(x=sapply(steps, c), y=sapply(values,c), min=sapply(min, c), max=sapply(max, c), median=sapply(median, c))
  df_final
}

subplot_variable_canopy <- function(var_name, start, end, liste_df_ent_1, liste_df_ent_2, liste_legend_name, plottitle, ytitle, lai_surfsol=1)
{
  df_1 <- plot_dataframe_canopy(var_name, liste_df_ent_1, liste_legend_name, start, end)
  df_2 <- plot_dataframe_canopy(var_name, liste_df_ent_2, liste_legend_name, start, end)

  df <- df_1
  df$values <- df$values + df_2$values

  if (lai_surfsol != 1)
  {
    df$values <- sapply(df$values, function(x){x/lai_surfsol})
  }

  p <- ggplot(data=df, aes(x=x, y=values, group=light)) +
    geom_line(aes(color=light)) +
    ggtitle(plottitle) +
    xlab("Jour de l'année") +
    ylab(ytitle) +
    labs(color="Modèle")
  p <- p+theme(plot.title = element_text(size=11))
  p
}

plot_variable_canopy <- function(var_name, start, end, listes_situation_1, listes_situation_2, liste_legend_name, ytitle, globaltitle, passive=FALSE, lai_surfsol=1)
{
  if (passive)
  {
    plottitle <- 'Photomorpho activée, Ramifications activées, RATP passif'
    p1 <- subplot_variable_canopy(var_name, start, end, listes_situation_1[[1]], listes_situation_2[[1]], liste_legend_name, plottitle, ytitle, lai_surfsol)
  }


  plottitle <- 'Photomorpho activée, Ramifications désactivées, RATP actif'
  p2 <- subplot_variable_canopy(var_name, start, end, listes_situation_1[[2]], listes_situation_2[[2]], liste_legend_name, plottitle, ytitle, lai_surfsol)

  plottitle <- 'Photomorpho activée, Ramifications activées, RATP actif'
  p3 <- subplot_variable_canopy(var_name, start, end, listes_situation_1[[3]], listes_situation_2[[3]], liste_legend_name, plottitle, ytitle, lai_surfsol)

  plottitle <- 'Photomorpho désactivée, RATP actif'
  p4 <- subplot_variable_canopy(var_name, start, end, listes_situation_1[[4]], listes_situation_2[[4]], liste_legend_name, plottitle, ytitle, lai_surfsol)

  ## agencement des graphes
  if (passive)
  {
    p_init <- plot_grid(
      p1 + theme(legend.position = "none"),
      p2 + theme(legend.position = "none"),
      p3 + theme(legend.position = "none"),
      p4 + theme(legend.position = "none"),
      ncol=2, nrow=2, labels=c("A","B","C","D")
    )
  }
  else
  {
    p_init <- plot_grid(
      p2 + theme(legend.position = "none"),
      p3 + theme(legend.position = "none"),
      p4 + theme(legend.position = "none"),
      ncol=2, nrow=2, labels=c("A","B","C")
    )
  }


  legend <- get_legend(
    p2 + theme(legend.box.margin = margin(0,0,0,12))
  )

  p_legend <- plot_grid(p_init, legend, nrow=1, rel_widths = c(5, 1))


  title <- ggdraw() + draw_label(globaltitle, fontface='bold')

  p_final <- plot_grid(title, p_legend,ncol=1, rel_heights = c(0.1, 1))
  p_final
}

subplot_variable_entity <- function(var_name, start, end, liste_df_ent, liste_legend_name, plottitle, ytitle, lai_surfsol=1)
{
  df <- plot_dataframe_canopy(var_name, liste_df_ent, liste_legend_name, start, end)

  if (lai_surfsol != 1)
  {
    df$values <- sapply(df$values, function(x){x/lai_surfsol})
  }

  p <- ggplot(data=df, aes(x=x, y=values, group=light)) +
    geom_line(aes(color=light)) +
    ggtitle(plottitle) +
    xlab("Jour de l'année") +
    ylab(ytitle) +
    labs(color="Modèle")
  p <- p+theme(plot.title = element_text(size=11))
  p
}

plot_variable_entity <- function(var_name, start, end, listes_situation, liste_legend_name, ytitle, globaltitle, passive=FALSE, lai_surfsol=1)
{
  if (passive)
  {
    plottitle <- 'Photomorpho activée, Ramifications activées, RATP passif'
    p1 <- subplot_variable_entity(var_name, start, end, listes_situation[[1]], liste_legend_name, plottitle, ytitle, lai_surfsol)
  }


  plottitle <- 'Photomorpho activée, Ramifications désactivées, RATP actif'
  p2 <- subplot_variable_entity(var_name, start, end, listes_situation[[2]], liste_legend_name, plottitle, ytitle, lai_surfsol)

  plottitle <- 'Photomorpho activée, Ramifications activées, RATP actif'
  p3 <- subplot_variable_entity(var_name, start, end, listes_situation[[3]], liste_legend_name, plottitle, ytitle, lai_surfsol)

  plottitle <- 'Photomorpho désactivée, RATP actif'
  p4 <- subplot_variable_entity(var_name, start, end, listes_situation[[4]], liste_legend_name, plottitle, ytitle, lai_surfsol)

  ## agencement des graphes
  if (passive)
  {
    p_init <- plot_grid(
      p1 + theme(legend.position = "none"),
      p2 + theme(legend.position = "none"),
      p3 + theme(legend.position = "none"),
      p4 + theme(legend.position = "none"),
      ncol=2, nrow=2, labels=c("A","B","C","D")
    )
  }
  else
  {
    p_init <- plot_grid(
      p2 + theme(legend.position = "none"),
      p3 + theme(legend.position = "none"),
      p4 + theme(legend.position = "none"),
      ncol=2, nrow=2, labels=c("A","B","C")
    )
  }


  legend <- get_legend(
    p2 + theme(legend.box.margin = margin(0,0,0,12))
  )

  p_legend <- plot_grid(p_init, legend, nrow=1, rel_widths = c(5, 1))


  title <- ggdraw() + draw_label(globaltitle, fontface='bold')

  p_final <- plot_grid(title, p_legend,ncol=1, rel_heights = c(0.1, 1))
  p_final
}

diff_dataframe_toto_legume <- function(df1, df2, passive=FALSE)
{
  df_final <- df1
  if (passive)
  {
    for (i in 3:ncol(df2))
    {
      df_final[df_final$V1 == "epsi" ,][,i] <- as.numeric(as.character(df_final[df_final$V1 == "epsi" ,][,i])) - as.numeric(as.character(df2[df2$V1 == "epsi" ,][,i]))
      df_final[df_final$V1 == "PARiPlante" ,][,i] <- as.numeric(as.character(df_final[df_final$V1 == "PARiPlante" ,][,i])) - as.numeric(as.character(df2[df2$V1 == "PARiPlante" ,][,i]))
    }
  }
  else
  {
    for (i in 3:ncol(df2))
    {
      df_final[2:nrow(df_final),i] <- as.numeric(as.character(df_final[2:nrow(df_final),i])) - as.numeric(as.character(df2[2:nrow(df_final),i]))
    }
  }

  df_final
}

subplot_variable_plant <- function(var_name, start, end, liste_df_ent, liste_legend_name, plantlist, plottitle, ytitle, lai_surfsol=1)
{
  df <- plot_dataframe_plant(var_name, liste_df_ent, liste_legend_name, plantlist, t_start=start, t_end=end)

  if (lai_surfsol != 1)
  {
    df$values <- sapply(df$values, function(x){x/lai_surfsol})
  }

  # for (j in 1:length(plantlist))
  # {
  #
  #     n <- paste(liste_legend_name[[1]], as.character(plantlist[[j]] -3), sep=" plante ")
  #     print(n)
  #     print(max(as.numeric(df[df$light==n,2])))
  # }

  p <- ggplot(data=df, aes(x=x, y=values, group=light)) +
    geom_line(aes(color=light)) +
    ggtitle(plottitle) +
    xlab("Jour de l'année") +
    ylab(ytitle) +
    labs(color="Modèle")
  p <- p+theme(plot.title = element_text(size=11))
  p
}

plot_variable_plant <- function(var_name, start, end, listes_situation, liste_legend_name, plantlist, ytitle, globaltitle, passive=FALSE, lai_surfsol=1)
{
  if (passive)
  {
    plottitle <- 'Photomorpho activée, Ramifications activées, RATP passif'
    p1 <- subplot_variable_plant(var_name, start, end, listes_situation[[1]], liste_legend_name, plantlist, plottitle, ytitle, lai_surfsol)
  }


  plottitle <- 'Photomorpho activée, Ramifications désactivées, RATP actif'
  p2 <- subplot_variable_plant(var_name, start, end, listes_situation[[2]], liste_legend_name, plantlist, plottitle, ytitle, lai_surfsol)

  plottitle <- 'Photomorpho activée, Ramifications activées, RATP actif'
  p3 <- subplot_variable_plant(var_name, start, end, listes_situation[[3]], liste_legend_name, plantlist, plottitle, ytitle, lai_surfsol)

  plottitle <- 'Photomorpho désactivée, RATP actif'
  p4 <- subplot_variable_plant(var_name, start, end, listes_situation[[4]], liste_legend_name, plantlist, plottitle, ytitle, lai_surfsol)

  ## agencement des graphes
  if (passive)
  {
    p_init <- plot_grid(
      p1 + theme(legend.position = "none"),
      p2 + theme(legend.position = "none"),
      p3 + theme(legend.position = "none"),
      p4 + theme(legend.position = "none"),
      ncol=2, nrow=2, labels=c("A","B","C","D")
    )
  }
  else
  {
    p_init <- plot_grid(
      p2 + theme(legend.position = "none"),
      p3 + theme(legend.position = "none"),
      p4 + theme(legend.position = "none"),
      ncol=2, nrow=2, labels=c("A","B","C")
    )
  }


  legend <- get_legend(
    p2 + theme(legend.box.margin = margin(0,0,0,12))
  )

  p_legend <- plot_grid(p_init, legend, nrow=1, rel_widths = c(5, 1))


  title <- ggdraw() + draw_label(globaltitle, fontface='bold')

  p_final <- plot_grid(title, p_legend,ncol=1, rel_heights = c(0.1, 1))
  p_final
}

max_per_plant <- function(var_name, df)
{
  plante <- list()
  maxvar <- list()
  for (i in 3:ncol(df))
  {
    maxvar <- append(maxvar, max(as.numeric(df[df$V1 == varname ,i])))
    plante <- append(plante, i)
  }


  df_final <- data.frame(plante=sapply(plante,c), max=sapply(maxvar,c))
  df_final
}

df_correlation_plante <- function(var_name, default_list, ratp_list, situation_list)
{
  lightmodel <- list()
  default <- list()
  situation <- list()
  for (i in 1 : length(default_list))
  {
    for (j in 3 : ncol(default_list[[i]]))
    {
      default <- append(default, as.numeric(default_list[[i]][default_list[[i]]$V1 == var_name,][,j]))
      lightmodel <- append(lightmodel, as.numeric(ratp_list[[i]][ratp_list[[i]]$V1 == var_name,][,j]))
      situation <- append(situation, rep(situation_list[[i]], length(ratp_list[[i]][ratp_list[[i]]$V1 == var_name,][,j])))
    }
  }

  df <- data.frame(default = sapply(default, c), lightmodel=sapply(lightmodel, c), situation=sapply(situation, c))
  df
}

df_correlation_plante_onestep <- function(var_name1, var_name2, default_list, model_list, step)
{
  lightmodel <- list()
  default <- list()
  varcolor <- list()
  espece <- list()
  for (i in 1 : length(default_list))
  {
    for (j in 3 : ncol(default_list[[i]]))
    {
      default <- append(default, as.numeric(default_list[[i]][default_list[[i]]$V1 == var_name1 & default_list[[i]]$V2 == step,][,j]))
      lightmodel <- append(lightmodel, as.numeric(model_list[[i]][model_list[[i]]$V1 == var_name1 & model_list[[i]]$V2 == step,][,j]))
      varcolor <- append(varcolor, as.numeric(default_list[[i]][default_list[[i]]$V1 == var_name2 & default_list[[i]]$V2 == step,][,j]))
      espece <- append(espece, rep(as.character(i), length(default_list[[i]][default_list[[i]]$V1 == var_name1 & default_list[[i]]$V2 == step,][,j])))
    }
  }

  df <- data.frame(default = sapply(default, c), lightmodel=sapply(lightmodel, c), var2=sapply(varcolor, c), espece=sapply(espece, c))
  df
}

plante_onestep_correlation <- function(varname1,
                                       varname2,
                                       step,
                                       default_list_case1, model_list_case1,
                                       default_list_case2, model_list_case2,
                                       xymax,
                                       legendname,
                                       case1name,
                                       case2name,
                                       surfsol=1)
{
  df_xy <- data.frame(x=c(0,xymax),y=c(0,xymax))

  df <- df_correlation_plante_onestep(varname1, varname2, default_list_case1, model_list_case1, step)
  df$var2 <- sapply(df$var2, function(x){x/surfsol})

  p1 <- ggplot() +
    geom_point(data=df, aes(x=default, y=lightmodel, color=var2, shape=espece), size = 3) +
    geom_line(data=df_xy,aes(x=x,y=y)) +
    scale_colour_gradientn(colours=rainbow(2)) +
    labs(x="Défaut l-egume",y="CARIBU Actif",title=case1name, color=legendname)
  p1

  df <- df_correlation_plante_onestep(varname1, varname2, default_list_case2, model_list_case2, step)
  df$var2 <- sapply(df$var2, function(x){x/surfsol})

  p2 <- ggplot() +
    geom_point(data=df, aes(x=default, y=lightmodel, color=var2, shape=espece), size = 3) +
    geom_line(data=df_xy,aes(x=x,y=y)) +
    scale_colour_gradientn(colours=rainbow(2)) +
    labs(x="Défaut l-egume",y="CARIBU Actif",title=case2name, color=legendname)
  p2

  p_init <- plot_grid(
    p1,
    p2,
    ncol=2, nrow=1, labels=c("A","B")
  )

  globaltitle <- paste("Pas de temps :", as.character(step), sep=" ")
  title <- ggdraw() + draw_label(globaltitle, fontface='bold')

  p_final <- plot_grid(title, p_init,ncol=1, rel_heights = c(0.1, 1))
  p_final
}

