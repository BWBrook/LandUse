### Food demand analysis ###
rm(list=ls(all=TRUE))
source('C:/GitHub/LandUseChange/AIC_BIC_functions.r', local=T)
source('C:/GitHub/LandUseChange/ArtOfRProg.r', local=T)

setwd("C:/GitHub/LandUseChange")

dat.kcal <- read.table("FAOkcal.csv", header=T, sep=",", row.names=1)
dat.gdp <- read.table("FAOgdp.csv", header=T, sep=",", row.names=1)

fix(dat.kcal) ## Check imported file integrity
fix(dat.gdp)

plot(dat.gdp$Y1961, dat.kcal$Y1961)
plot(dat.gdp$Y2008, dat.kcal$Y2008)
plot(as.numeric(dat.gdp["Australia",]), as.numeric(dat.kcal["Australia",]))

cor(as.numeric(dat.gdp["Australia",]), as.numeric(dat.kcal["Australia",]), method="kendall")
udcorr(as.numeric(dat.gdp["Australia",]), as.numeric(dat.kcal["Australia",]))

