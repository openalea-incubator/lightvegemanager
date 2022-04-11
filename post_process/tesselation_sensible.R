library(ggplot2)

# se place dans le répertoire
setwd("C:/Users/mwoussen/cdd/codes/vegecouplelight/outputs/dynamic_cn_wheat_csv")


filenames <- paste(paste("tesselation_sensibility_",1,sep=""), "csv", sep=".")
for (i in 2:12)
{
  filenames <- append(filenames, paste(paste("tesselation_sensibility_",i,sep=""), "csv", sep="."))
}

data <- lapply(filenames, read.csv)


# boite à moustache pour chaque x
# x : niveau tesselation
# y : différence (PARi_CARIBU - PARa_RATP)_t 
df_jour <- data[[1]][data[[1]]["PAR.input"] > 0 ,]
diff <- df_jour["PARi.CARIBU"] - df_jour["PARa.RATP"]
for (i in 2:12)
{
  df_jour <- data[[i]][data[[i]]["PAR.input"] > 0 ,]
  diff <- append(diff, df_jour["PARi.CARIBU"] - df_jour["PARa.RATP"])
}

# refonte des données en dataframe
# "values" "level"
df_diff <- data.frame(val=diff[[1]], level=list(rep(1, length(diff[[1]])))[[1]])
for(i in 2:12)
{
  df_diff <- rbind(df_diff, data.frame(val=diff[[i]], level=list(rep(i, length(diff[[i]])))[[1]]))
}

df_diff$level <- as.factor(df_diff$level)

# plot 
p <- ggplot(df_diff, aes(x=level, y=val)) + geom_boxplot() + 
        labs(title="Comparaison sur 3 jours", x="Niveau de tesselation des triangles", y="différence PAR incident CARIBU - PARa RATP (µmol.m-2.s-1)")
p