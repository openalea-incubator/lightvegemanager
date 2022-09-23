library(ggplot2)

# se place dans le dossier
setwd('C:/Users/mngauthier/mwoussen_2022/postprocess_legume_ratp')

legume_default <- read.table('toto_17112_l-egume_seul.csv', sep=';',stringsAsFactors = FALSE)
legume_ratp <- read.table('toto_17112_l-egume_ratp_1.csv', sep=';', stringsAsFactors = FALSE)
legume_ratp_passive <- read.table('outputs_ratp_passive.csv', sep=',', stringsAsFactors = TRUE)

# trace des graphes
# inputs : 
#   - liste de dataframe
#   - nom variable à comparer
#   - préciser si temps réduit
plot_dataframe_canopy <- function(variable, dataframe_list, name_list, t_start=-1, t_end=-1){
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
  df = data.frame(x=sapply(steps,c), values=sapply(values,c), light=sapply(names,c))
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

## Epsi ##
# df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp_passive), list("default", "ratp"))
df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp_passive), list("default", "ratp"))
# df <- plot_dataframe_canopy("epsi", list(legume_ratp_passive), list("ratp"))
df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp, legume_ratp_passive), list("default", "ratp active", "ratp passive"))

p <- ggplot(data=df, aes(x=x, y=values, group=light)) + 
      geom_line(aes(color=light)) +
      ggtitle('Epsi sommé sur le couvert') +
      xlab("DOY") +
      ylab("epsi")
p

## SurfPlante ##
# df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp_passive), list("default", "ratp"))
df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp_passive), list("default", "ratp"))
# df <- plot_dataframe_canopy("epsi", list(legume_ratp_passive), list("ratp"))
df <- plot_dataframe_canopy("epsi", list(legume_default, legume_ratp, legume_ratp_passive), list("default", "ratp active", "ratp passive"))

p <- ggplot(data=df, aes(x=x, y=values, group=light)) + 
  geom_line(aes(color=light)) +
  ggtitle('Epsi sommé sur le couvert') +
  xlab("DOY") +
  ylab("epsi")
p



      
