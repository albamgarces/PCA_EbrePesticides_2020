##Environment preparation
# setwd(".")
# dir.create("./PesticidesEbre_PCA")
# setwd("./PesticidesEbre_PCA")
# dir.create("data")
# dir.create("results")
# dir.create("tables&figures")

##Data preparation
library(readr)
pesticides <- as.data.frame(read_csv("data/pesticides_Ebro_2020.csv"))
str(pesticides)
pesticides$bay <- as.factor(pesticides$bay)
pesticides$channel <- as.factor(pesticides$channel)
names(pesticides)[1] <- "site"
row.names(pesticides)<- pesticides$site

##data exploration
str(pesticides)
summary(pesticides)

##descriptive statistic
#-< descargar statistics_pesticides.csv de github a la carpeta data
statistics <- as.data.frame(read_csv("data/statistics_pesticides.csv",
                                     na="empty"))
##pairwise correlations
#install.packages("corrplot")
source("http://www.sthda.com/upload/rquery_cormat.r")
require(stats)
pest_corr <- cor(pesticides[,4:38])
write.csv(pest_corr, "./results/corr_pesticides.csv")

det(pest_corr)
#-< near to 0, there're correlation

#ordered correlation plot
corr_plot <- rquery.cormat(pesticides[,4:38])

# png("./tables&figures/figure1.png")
# rquery.cormat(pesticides[,4:38])$sym
# dev.off()

#network plot correlations
library(qgraph)
pest_groups <- qgraph(pest_corr, graph="cor", layout="spring", threshold=0.6, title="Network plot correlations")

# png("./tables&figures/figure2.png")
# qgraph(pest_corr, graph="cor", layout="spring", threshold=0.6)
# dev.off()

#a ranked cross-correlations package
#devtools::install_github("laresbernardo/lares")
library(lares)
corr_cross(pesticides[,4:38],
           max_pvalue = 0.05,
           top=10)

# png("./tables&figures/figure3.png")
# corr_cross(pesticides[,4:38],
#            max_pvalue = 0.05,
#            top=10)
#dev.off()

##componentes principales
#calculamos para cada pesticida la varianza
round(apply(pesticides[,4:38], 2, var),3)
#Las dimensiones de las muestras varían mucho entre sí.
#escalamos los datos para estandarizar variables
acp <- prcomp(pesticides[,4:38],
              center=TRUE, #se resta la media de cada variable
              scale=TRUE) #divide por desv tipica
summary(acp)
plot(acp, type="l")

#barplots four firsts PC
par(mfrow=c(2,2))
barplot(acp$rotation[,1], las=2, cex.names = 0.4, main = "a) PC1: 22.1%")
barplot(acp$rotation[,2], las=2, cex.names = 0.4, main = "b) PC2: 15.1%")
barplot(acp$rotation[,3], las=2, cex.names = 0.4, main = "c) PC3: 12%")
barplot(acp$rotation[,4], las=2, cex.names = 0.4, main = "d) PC4: 9.7%")

# png("./tables&figures/figure4.png")
# par(mfrow=c(2,2))
# barplot(acp$rotation[,1], las=2, cex.names = 0.4, main = "a) PC1: 22.1%")
# barplot(acp$rotation[,2], las=2, cex.names = 0.4, main = "b) PC2: 15.1%")
# barplot(acp$rotation[,3], las=2, cex.names = 0.4, main = "c) PC3: 12%")
# barplot(acp$rotation[,4], las=2, cex.names = 0.4, main = "d) PC4: 9.7%")
# dev.off()

#general PCA biplot
library(factoextra)

fviz_pca_var(acp, axes=c(1,2), col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE,
             title="a")
fviz_pca_var(acp, axes=c(3,4), col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE,
             title="b")
# png("./tables&figures/figure5a.png")
# fviz_pca_var(acp, axes=c(1,2), col.var = "contrib", 
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#              repel = TRUE,
#              title="a")
# dev.off()
# png("./tables&figures/figure5b.png")
# fviz_pca_var(acp, axes=c(3,4), col.var = "contrib", 
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#              repel = TRUE,
#              title="b")
# dev.off()

#scores and loadings plot PC1 and PC2

bay <- as.numeric(pesticides$bay)
channel <- as.numeric(pesticides$channel)

plot(acp$x[,1],
     pch=c(16,17)[bay],
     col=c("green", "orange")[channel],
     ylab="PC1",
     main="a")
text(acp$x[,1],labels=pesticides$site, pos = c(3, 4)[channel], cex=0.5)
abline(h=0, v=0, lty=4, col="gray")

plot(acp$x[,2],
     pch=c(16,17)[bay],
     col=c("green", "orange")[channel],
     ylab="PC2",
     main="b")
text(acp$x[,2],labels=pesticides$site, pos = c(3, 4)[channel], cex=0.5)
abline(h=0, v=0, lty=4, col="gray")

# png("./tables&figures/figure6a.png")
# plot(acp$x[,1],
#      pch=c(16,17)[bay],
#      col=c("green", "orange")[channel],
#      ylab="PC1",
#      main="a")
# text(acp$x[,1],labels=pesticides$site, pos = c(3, 4)[channel], cex=0.5)
# abline(h=0, v=0, lty=4, col="gray")
# dev.off()
# png("./tables&figures/figure6b.png")
# plot(acp$x[,2],
#      pch=c(16,17)[bay],
#      col=c("green", "orange")[channel],
#      ylab="PC2",
#      main="b")
# text(acp$x[,2],labels=pesticides$site, pos = c(3, 4)[channel], cex=0.5)
# abline(h=0, v=0, lty=4, col="gray")
# dev.off()

#PC1 and PC2 bliplot. 
#Acidics, amilides, organophosphates, triazines
#removed assymetry quoficient negative

fviz_pca_biplot(acp, col.var = "contrib", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                select.var = list(name=c("2,4-D", "Bentazone", "MCPA", 
                                         "Diflufenican", "Propanil", 
                                         "Atrazine", "Cybutrine", "Simazine", "Terbutryn",
                                         "Azinphos ethyl", "Chlorfenviphos", "Diazinon", "Malaoxon")),
                repel = TRUE)

# png("./tables&figures/figure7a.png")
# fviz_pca_biplot(acp, col.var = "contrib",
#                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                 select.var = list(name=c("2,4-D", "Bentazone", "MCPA",
#                                          "Diflufenican", "Propanil",
#                                          "Atrazine", "Cybutrine", "Simazine", "Terbutryn",
#                                          "Azinphos ethyl", "Chlorfenviphos", "Diazinon", "Malaoxon")),
#                 repel = TRUE)
# dev.off()


#Acidics, amilides, organophosphates, triazines
#all
fviz_pca_biplot(acp, col.var = "contrib", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                select.var = list(name=c("2,4-D", "Bentazone", "MCPA", 
                                         "Diflufenican", "Propanil", 
                                         "Atrazine", "Cybutrine", "Simazine", "Terbuthylazine", "Terbutryn",
                                         "Azinphos ethyl", "Chlorfenviphos", "Chlorpyrifos", "Diazinon", "Malaoxon")),
                repel = TRUE)


#PC3 and PC4 bliplot. 
#Acidics, amilides, organophosphates, triazines
#removed assymetry quoficient negative
fviz_pca_biplot(acp, axes=c(3,4),
                col.var = "contrib", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                select.var = list(name=c("2,4-D", "Bentazone", "MCPA", 
                                         "Diflufenican", "Propanil", 
                                         "Atrazine", "Cybutrine", "Simazine", "Terbutryn",
                                         "Azinphos ethyl", "Chlorfenviphos", "Diazinon", "Malaoxon")),
                repel = TRUE)

# png("./tables&figures/figure7b.png")
# fviz_pca_biplot(acp, axes=c(3,4),
#                 col.var = "contrib", 
#                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                 select.var = list(name=c("2,4-D", "Bentazone", "MCPA", 
#                                          "Diflufenican", "Propanil", 
#                                          "Atrazine", "Cybutrine", "Simazine", "Terbutryn",
#                                          "Azinphos ethyl", "Chlorfenviphos", "Diazinon", "Malaoxon")),
#                 repel = TRUE)
# dev.off()

#Acidics, amilides, organophosphates, triazines
#all
fviz_pca_biplot(acp, axes=c(3,4),
                col.var = "contrib", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                select.var = list(name=c("2,4-D", "Bentazone", "MCPA", 
                                         "Diflufenican", "Propanil", 
                                         "Atrazine", "Cybutrine", "Simazine", "Terbuthylazine", "Terbutryn",
                                         "Azinphos ethyl", "Chlorfenviphos", "Chlorpyrifos", "Diazinon", "Malaoxon")),
                repel = TRUE)