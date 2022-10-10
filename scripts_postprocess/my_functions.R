### TRAITEMENT SORTIES L-EGUME ###

# courbes somme sur le couvert
# inputs : 
#   - liste de dataframe
#   - nom variable Ã  comparer
#   - prÃ©ciser si temps rÃ©duit
plot_dataframe_canopy <- function(variable, dataframe_list, name_list, t_start=-1, t_end=-1)
{
  # nouvelle liste colonne : steps, value, name
  values <- list()
  steps <- list()
  names <- list()
  for (i in 1:length(dataframe_list))
  {
    if (ncol(dataframe_list[[i]]) > 5)
    {
      values <- append(values, rowSums(sapply(dataframe_list[[i]][dataframe_list[[i]]$V1 == variable ,][,3:34], function(x) as.numeric(as.character(x)))))
    }else
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
      names <- append(names, rep(paste(name_list[[i]], as.character(plant_list[[j]]), sep=" plante "), 
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
      steps <- append(steps, as.character(i))
      names <- append(names, func_names[[k]])
      
    }
  }
  
  df_final <- data.frame(x=sapply(steps,c), values=sapply(values,c), names=sapply(names, c))
  df_final
}
