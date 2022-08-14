# library(R2jags)
# library(coda)
# library(boot)
library(INLA)
library(here)
library(raster)
library(boot)
library(rgdal)
library(MetBrewer)
library(GISTools)
##2001-2010

rm(list=ls()) ## <----- add this line to clean up the workspace and make sure old objects don't accidentally carry over

#set_here("~/Box/GM - Dispar Moth/Jeffress Spread Forecasting") #<------ this line isn't necessary
## NOTE: some scripts are broken because data paths point all the way back to the root
## on a specific machine. Please use here() and relative paths as below.
here()


habitat<-raster(here("Data/forest_canopy_5km.tif"))
##connect 2001-2010
connectivity01<-raster(here("JuliaScripts/Full Extent/2001_5km_output/cum_currmap.tif"))
connectivity02<-raster(here("JuliaScripts/Full Extent/2002_5km_output/cum_currmap.tif"))
connectivity03<-raster(here("JuliaScripts/Full Extent/2003_5km_output/cum_currmap.tif"))
connectivity04<-raster(here("JuliaScripts/Full Extent/2004_5km_output/cum_currmap.tif"))
connectivity05<-raster(here("JuliaScripts/Full Extent/2005_5km_output/cum_currmap.tif"))
connectivity06<-raster(here("JuliaScripts/Full Extent/2006_5km_output/cum_currmap.tif"))
connectivity07<-raster(here("JuliaScripts/Full Extent/2007_5km_output/cum_currmap.tif"))
connectivity08<-raster(here("JuliaScripts/Full Extent/2008_5km_output/cum_currmap.tif"))
connectivity09<-raster(here("JuliaScripts/Full Extent/2009_5km_output/cum_currmap.tif"))
connectivity10<-raster(here("JuliaScripts/Full Extent/2010_5km_output/cum_currmap.tif"))
connectivity11<-raster(here("JuliaScripts/Full Extent/2011_5km_output/cum_currmap.tif"))
connectivity12<-raster(here("JuliaScripts/Full Extent/2012_5km_output/cum_currmap.tif"))
connectivity13<-raster(here("JuliaScripts/Full Extent/2013_5km_output/cum_currmap.tif"))
connectivity14<-raster(here("JuliaScripts/Full Extent/2014_5km_output/cum_currmap.tif"))
connectivity15<-raster(here("JuliaScripts/Full Extent/2015_5km_output/cum_currmap.tif"))
connectivity16<-raster(here("JuliaScripts/Full Extent/2016_5km_output/cum_currmap.tif"))
connectivity17<-raster(here("JuliaScripts/Full Extent/2017_5km_output/cum_currmap.tif"))
connectivity18<-raster(here("JuliaScripts/Full Extent/2018_5km_output/cum_currmap.tif"))
connectivity19<-raster(here("JuliaScripts/Full Extent/2019_5km_output/cum_currmap.tif"))
##theta 2002-2011
theta02<-raster(here("Data/Full Extent/GridTraps2002_full_presence3.tif"))
theta03<-raster(here("Data/Full Extent/GridTraps2003_full_presence3.tif"))
theta04<-raster(here("Data/Full Extent/GridTraps2004_full_presence3.tif"))
theta05<-raster(here("Data/Full Extent/GridTraps2005_full_presence3.tif"))
theta06<-raster(here("Data/Full Extent/GridTraps2006_full_presence3.tif"))
theta07<-raster(here("Data/Full Extent/GridTraps2007_full_presence3.tif"))
theta08<-raster(here("Data/Full Extent/GridTraps2008_full_presence3.tif"))
theta09<-raster(here("Data/Full Extent/GridTraps2009_full_presence3.tif"))
theta10<-raster(here("Data/Full Extent/GridTraps2010_full_presence3.tif"))
theta11<-raster(here("Data/Full Extent/GridTraps2011_full_presence3.tif"))
theta12<-raster(here("Data/Full Extent/GridTraps2012_full_presence3.tif"))
theta13<-raster(here("Data/Full Extent/GridTraps2013_full_presence3.tif"))
theta14<-raster(here("Data/Full Extent/GridTraps2014_full_presence3.tif"))
theta15<-raster(here("Data/Full Extent/GridTraps2015_full_presence3.tif"))
theta16<-raster(here("Data/Full Extent/GridTraps2016_full_presence3.tif"))
theta17<-raster(here("Data/Full Extent/GridTraps2017_full_presence3.tif"))
theta18<-raster(here("Data/Full Extent/GridTraps2018_full_presence3.tif"))
theta19<-raster(here("Data/Full Extent/GridTraps2019_full_presence3.tif"))
nn <- length(values(habitat)) * 10 #x2 for two years


### format as a data frame

ten.df <- data.frame(theta = c(values(theta02)
                               ,values(theta03)
                               ,values(theta04)
                               ,values(theta05)
                               ,values(theta06)
                               ,values(theta07)
                               ,values(theta08)
                               ,values(theta09)
                               ,values(theta10)
                               ,values(theta11)
),
habitat = rep(values(habitat), times=10),
connect = log10(c(values(connectivity01) ## <--- NOTE addition of log10 here
                  ,values(connectivity02)
                  ,values(connectivity03)
                  ,values(connectivity04)
                  ,values(connectivity05)
                  ,values(connectivity06)
                  ,values(connectivity07)
                  ,values(connectivity08)
                  ,values(connectivity09)
                  ,values(connectivity10))
                +1)

)
ten.df <- ten.df[!is.na(ten.df$theta),]
ten.df <- ten.df[complete.cases(ten.df),]

### THIS SECTION IS NEW AND DIFFERENT FROM JAGS


## write down the model formula

formula = theta ~ connect + habitat

result = inla(formula, family="binomial", Ntrials=1, data=ten.df,
              control.family=list(link='logit'), verbose=TRUE,
              control.compute = list(waic=T))

summary(result)


## prediction time!

## logit(psi[t]) = b_0 + b_connect*connect[t] + b_habitat*habitat


predictions <- function(result, connect, habitat, psi.crit=0.5){ ##<------ I fixed an error in this line
  #params  =  a named list of parameter values for b_0, b_connect, etc.
  #connect = a vector of connectivity values
  #habitat = a vector of habitat values
  b_0 <- result$summary.fixed[1,1] ## <------- the next three lines are revised to use the output format from INLA
  b_connect <- result$summary.fixed[2,1]
  b_habitat <- result$summary.fixed[3,1]
  
  psi <- inv.logit(b_0 + b_connect*connect + b_habitat*habitat)
  theta <- ifelse(psi >= psi.crit, 1, 0)
  
  return(list(psi=psi, theta=theta))
}


## NOTE new log10 on connectivity values
### prediction accuracy for the year 2011
preds1 <- predictions(result, connect=log10(values(connectivity10)+1), habitat=values(habitat))
## true positive rate
true2011 <- sum(preds1$theta==1 & values(theta11)==1, na.rm=TRUE)/sum(values(theta11)==1, na.rm=TRUE)
true2011
## false positive rate
falsepos2011 <- sum(preds1$theta==1 & values(theta11)==0, na.rm=TRUE)/sum(values(theta11)==0, na.rm=TRUE)
falsepos2011
## total accuracy
total2011 <- mean(preds1$theta==values(theta11), na.rm=TRUE)
total2011

library(prettymapr)
## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta11))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2011, proj4string(theta11))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds1$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds1$theta==1 & values(theta11)==1] <- 1
vectorAll[preds1$theta==1 & values(theta11)==0] <- 2
vectorAll[preds1$theta==0 & values(theta11)==1] <- 3
vectorAll[preds1$theta==0 & values(theta11)==0] <- 4
## create a raster of said accuracy types
raster2011<-raster(matrix(vectorAll, nrow(theta11), ncol(theta11), byrow=T), crs=proj4string(theta11))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2011) <- e
raster2011_ext <- setExtent(raster2011, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2011_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2011 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)

### prediction accuracy for the year 2012
preds2 <- predictions(result, connect=log10(values(connectivity11)+1), habitat=values(habitat))
## true positive rate
true2012 <- sum(preds2$theta==1 & values(theta12)==1, na.rm=TRUE)/sum(values(theta12)==1, na.rm=TRUE)
true2012
## false positive rate
falsepos2012 <- sum(preds2$theta==1 & values(theta12)==0, na.rm=TRUE)/sum(values(theta12)==0, na.rm=TRUE)
falsepos2012
## total accuracy
total2012 <- mean(preds2$theta==values(theta12), na.rm=TRUE)
total2012

## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta12))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2012, proj4string(theta12))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds2$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds2$theta==1 & values(theta12)==1] <- 1
vectorAll[preds2$theta==1 & values(theta12)==0] <- 2
vectorAll[preds2$theta==0 & values(theta12)==1] <- 3
vectorAll[preds2$theta==0 & values(theta12)==0] <- 4
## create a raster of said accuracy types
raster2012<-raster(matrix(vectorAll, nrow(theta12), ncol(theta12), byrow=T), crs=proj4string(theta12))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2012) <- e
raster2012_ext <- setExtent(raster2012, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2012_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2012 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)

### prediction accuracy for the year 2013
preds3 <- predictions(result, connect=log10(values(connectivity12)+1), habitat=values(habitat))
## true positive rate
true2013 <- sum(preds3$theta==1 & values(theta13)==1, na.rm=TRUE)/sum(values(theta13)==1, na.rm=TRUE)
true2013
## false positive rate
falsepos2013 <- sum(preds3$theta==1 & values(theta13)==0, na.rm=TRUE)/sum(values(theta13)==0, na.rm=TRUE)
falsepos2013
## total accuracy
total2013 <- mean(preds3$theta==values(theta13), na.rm=TRUE)
total2013

## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta13))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2013, proj4string(theta13))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds3$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds3$theta==1 & values(theta13)==1] <- 1
vectorAll[preds3$theta==1 & values(theta13)==0] <- 2
vectorAll[preds3$theta==0 & values(theta13)==1] <- 3
vectorAll[preds3$theta==0 & values(theta13)==0] <- 4
## create a raster of said accuracy types
raster2013<-raster(matrix(vectorAll, nrow(theta13), ncol(theta13), byrow=T), crs=proj4string(theta13))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2013) <- e
raster2013_ext <- setExtent(raster2013, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2013_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2013 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)

### prediction accuracy for the year 2014
preds4 <- predictions(result, connect=log10(values(connectivity13)+1), habitat=values(habitat))
## true positive rate
true2014 <- sum(preds4$theta==1 & values(theta14)==1, na.rm=TRUE)/sum(values(theta14)==1, na.rm=TRUE)
true2014
## false positive rate
falsepos2014 <- sum(preds4$theta==1 & values(theta14)==0, na.rm=TRUE)/sum(values(theta14)==0, na.rm=TRUE)
falsepos2014
## total accuracy
total2014 <- mean(preds4$theta==values(theta14), na.rm=TRUE)
total2014

## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta14))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2014, proj4string(theta14))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds4$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds4$theta==1 & values(theta14)==1] <- 1
vectorAll[preds4$theta==1 & values(theta14)==0] <- 2
vectorAll[preds4$theta==0 & values(theta14)==1] <- 3
vectorAll[preds4$theta==0 & values(theta14)==0] <- 4
## create a raster of said accuracy types
raster2014<-raster(matrix(vectorAll, nrow(theta14), ncol(theta14), byrow=T), crs=proj4string(theta14))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2014) <- e
raster2014_ext <- setExtent(raster2014, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2014_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2014 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)

### prediction accuracy for the year 2015
preds5 <- predictions(result, connect=log10(values(connectivity14)+1), habitat=values(habitat))
## true positive rate
true2015 <- sum(preds5$theta==1 & values(theta15)==1, na.rm=TRUE)/sum(values(theta15)==1, na.rm=TRUE)
true2015
## false positive rate
falsepos2015 <- sum(preds5$theta==1 & values(theta15)==0, na.rm=TRUE)/sum(values(theta15)==0, na.rm=TRUE)
falsepos2015
## total accuracy
total2015 <- mean(preds5$theta==values(theta15), na.rm=TRUE)
total2015

## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta15))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2015, proj4string(theta15))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds5$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds5$theta==1 & values(theta15)==1] <- 1
vectorAll[preds5$theta==1 & values(theta15)==0] <- 2
vectorAll[preds5$theta==0 & values(theta15)==1] <- 3
vectorAll[preds5$theta==0 & values(theta15)==0] <- 4
## create a raster of said accuracy types
raster2015<-raster(matrix(vectorAll, nrow(theta15), ncol(theta15), byrow=T), crs=proj4string(theta15))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2015) <- e
raster2015_ext <- setExtent(raster2015, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2015_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2015 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)

### prediction accuracy for the year 2016
preds6 <- predictions(result, connect=log10(values(connectivity15)+1), habitat=values(habitat))
## true positive rate
true2016 <- sum(preds6$theta==1 & values(theta16)==1, na.rm=TRUE)/sum(values(theta16)==1, na.rm=TRUE)
true2016
## false positive rate
falsepos2016 <- sum(preds6$theta==1 & values(theta16)==0, na.rm=TRUE)/sum(values(theta16)==0, na.rm=TRUE)
falsepos2016
## total accuracy
total2016 <- mean(preds6$theta==values(theta16), na.rm=TRUE)
total2016

## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta16))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2016, proj4string(theta16))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds6$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds6$theta==1 & values(theta16)==1] <- 1
vectorAll[preds6$theta==1 & values(theta16)==0] <- 2
vectorAll[preds6$theta==0 & values(theta16)==1] <- 3
vectorAll[preds6$theta==0 & values(theta16)==0] <- 4
## create a raster of said accuracy types
raster2016<-raster(matrix(vectorAll, nrow(theta16), ncol(theta16), byrow=T), crs=proj4string(theta16))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2016) <- e
raster2016_ext <- setExtent(raster2016, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2016_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2016 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)

### prediction accuracy for the year 2017
preds7 <- predictions(result, connect=log10(values(connectivity16)+1), habitat=values(habitat))
## true positive rate
true2017 <- sum(preds7$theta==1 & values(theta17)==1, na.rm=TRUE)/sum(values(theta17)==1, na.rm=TRUE)
true2017
## false positive rate
falsepos2017 <- sum(preds7$theta==1 & values(theta17)==0, na.rm=TRUE)/sum(values(theta17)==0, na.rm=TRUE)
falsepos2017
## total accuracy
total2017 <- mean(preds7$theta==values(theta17), na.rm=TRUE)
total2017

## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta17))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2017, proj4string(theta17))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds7$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds7$theta==1 & values(theta17)==1] <- 1
vectorAll[preds7$theta==1 & values(theta17)==0] <- 2
vectorAll[preds7$theta==0 & values(theta17)==1] <- 3
vectorAll[preds7$theta==0 & values(theta17)==0] <- 4
## create a raster of said accuracy types
raster2017<-raster(matrix(vectorAll, nrow(theta17), ncol(theta17), byrow=T), crs=proj4string(theta17))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2017) <- e
raster2017_ext <- setExtent(raster2017, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2017_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2017 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)

### prediction accuracy for the year 2018
preds8 <- predictions(result, connect=log10(values(connectivity17)+1), habitat=values(habitat))
## true positive rate
true2018 <- sum(preds8$theta==1 & values(theta18)==1, na.rm=TRUE)/sum(values(theta18)==1, na.rm=TRUE)
true2018
## false positive rate
falsepos2018 <- sum(preds8$theta==1 & values(theta18)==0, na.rm=TRUE)/sum(values(theta18)==0, na.rm=TRUE)
falsepos2018
## total accuracy
total2018 <- mean(preds8$theta==values(theta18), na.rm=TRUE)
total2018

## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta18))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2018, proj4string(theta18))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds8$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds8$theta==1 & values(theta18)==1] <- 1
vectorAll[preds8$theta==1 & values(theta18)==0] <- 2
vectorAll[preds8$theta==0 & values(theta18)==1] <- 3
vectorAll[preds8$theta==0 & values(theta18)==0] <- 4
## create a raster of said accuracy types
raster2018<-raster(matrix(vectorAll, nrow(theta18), ncol(theta18), byrow=T), crs=proj4string(theta18))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2018) <- e
raster2018_ext <- setExtent(raster2018, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2018_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2018 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)

### prediction accuracy for the year 2019
preds9 <- predictions(result, connect=log10(values(connectivity18)+1), habitat=values(habitat))
## true positive rate
true2019 <- sum(preds9$theta==1 & values(theta19)==1, na.rm=TRUE)/sum(values(theta19)==1, na.rm=TRUE)
true2019
## false positive rate
falsepos2019 <- sum(preds9$theta==1 & values(theta19)==0, na.rm=TRUE)/sum(values(theta19)==0, na.rm=TRUE)
falsepos2019
## total accuracy
total2019 <- mean(preds9$theta==values(theta19), na.rm=TRUE)
total2019

## read the states shape file
shp <- readOGR( "./Data/Raw/statesp020.shp") 
## transform the shape file to the projection of our traps
trans_shp <- spTransform(shp, proj4string(theta19))
## also transform the accuracy types to match this extent
trans_rast <- spTransform(raster2019, proj4string(theta19))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
## crop the shape file to the extent of the accuracy types raster
shp_crop <- crop(trans_shp, e)

## create an empty vector
vectorAll <- rep_len(NA, length(preds9$theta))
## fill the vector in with conditional codes for accuracy types
vectorAll[preds9$theta==1 & values(theta19)==1] <- 1
vectorAll[preds9$theta==1 & values(theta19)==0] <- 2
vectorAll[preds9$theta==0 & values(theta19)==1] <- 3
vectorAll[preds9$theta==0 & values(theta19)==0] <- 4
## create a raster of said accuracy types
raster2019<-raster(matrix(vectorAll, nrow(theta19), ncol(theta19), byrow=T), crs=proj4string(theta19))
e <- extent(-1e+05, 1900000, 1300000, 2950000)
extent(raster2019) <- e
raster2019_ext <- setExtent(raster2019, e, keepres=TRUE)
#plot raster with categorical legend!
plot(raster2019_ext,
     legend = FALSE,
     col = c("#225bb2","#973d21","#ee956a","#7db0ea"),
     axes = FALSE,
     par(mar=c(1,1,1,1))
)

legend("topright",
       title = "2019 Performance", 
       legend = c("True Positive", "False Positive",
                  "False Negative", "True Negative"),
       fill = c("#225bb2","#973d21","#ee956a","#7db0ea"),
       cex=1.10,
       bty = "n")
subset <- subset(shp_crop, STATE=="Virginia"| STATE=="Wisconsin"
                 | STATE=="Minnesota"| STATE=="Ohio"| STATE=="Indiana"| STATE=="North Carolina"
                 | STATE=="Tennessee"| STATE=="Kentucky"| STATE=="Illinois"| STATE=="West Virginia"
                 | STATE=="Iowa"| STATE=="Missouri")
plot(subset, add=TRUE)
addscalebar(widthhint=.30)
## plots #######################################################################
year <- c(2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)
accuracy <- c(total2011, total2012, total2013, total2014, total2015, total2016, total2017, total2018, total2019)
accuracypos <- c(true2011, true2012, true2013, true2014, true2015, true2016, true2017, true2018, true2019)
falsepos <- c(falsepos2011, falsepos2012, falsepos2013, falsepos2014, falsepos2015, falsepos2016, falsepos2017, falsepos2018, falsepos2019)
## simple scatterplot
attach(mtcars)
par(mfrow=c(1,3))
plot(year, accuracy, main="5km Radius, Presence 3", 
     xlab="Year ", ylab="Total Accuracy", pch=19)
plot(year, accuracypos, main="5km Radius, Presence 3", 
     xlab="Year ", ylab="True Positive", pch=19)
plot(year, falsepos, main="5km Radius, Presence 3", 
     xlab="Year ", ylab="False Positive", pch=19)
mean(accuracy)
mean(accuracypos)
mean(falsepos)
## NOTE: write these outputs to new directory ./Data/Model Predictions
## Write rasters to images
psi.pred.raster11<-raster(matrix(preds1$psi, nrow(theta11), ncol(theta11), byrow=T), crs=proj4string(theta11))
writeRaster(x=psi.pred.raster11, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2011_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

psi.pred.raster12<-raster(matrix(preds2$psi, nrow(theta12), ncol(theta12), byrow=T), crs=proj4string(theta12))
writeRaster(x=psi.pred.raster12, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2012_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

psi.pred.raster13<-raster(matrix(preds3$psi, nrow(theta13), ncol(theta13), byrow=T), crs=proj4string(theta13))
writeRaster(x=psi.pred.raster13, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2013_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

psi.pred.raster14<-raster(matrix(preds4$psi, nrow(theta14), ncol(theta14), byrow=T), crs=proj4string(theta14))
writeRaster(x=psi.pred.raster14, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2014_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

psi.pred.raster15<-raster(matrix(preds5$psi, nrow(theta15), ncol(theta15), byrow=T), crs=proj4string(theta15))
writeRaster(x=psi.pred.raster15, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2015_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

psi.pred.raster16<-raster(matrix(preds6$psi, nrow(theta16), ncol(theta16), byrow=T), crs=proj4string(theta16))
writeRaster(x=psi.pred.raster16, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2016_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

psi.pred.raster17<-raster(matrix(preds7$psi, nrow(theta17), ncol(theta17), byrow=T), crs=proj4string(theta17))
writeRaster(x=psi.pred.raster17, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2017_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

psi.pred.raster18<-raster(matrix(preds8$psi, nrow(theta18), ncol(theta18), byrow=T), crs=proj4string(theta18))
writeRaster(x=psi.pred.raster18, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2018_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

psi.pred.raster19<-raster(matrix(preds9$psi, nrow(theta19), ncol(theta19), byrow=T), crs=proj4string(theta19))
writeRaster(x=psi.pred.raster19, filename=here("./Data/Reduced/2001_2010_5km_5km_3_2019_INLA.tif"), 
            format="GTiff", overwrite=TRUE)

## Plots!!
par(mfrow = c(3, 2))
plot(psi.pred.raster11, main="5km, 3, predictions 2011")
plot(theta11, main="5km, 3, actual 2011")

plot(psi.pred.raster12, main="5km, 3, predictions 2012")
plot(theta12, main="5km, 3, actual 2012")

plot(psi.pred.raster13, main="5km, 3, predictions 2013")
plot(theta13, main="5km, 3, actual 2013")

plot(psi.pred.raster14, main="5km, 3, predictions 2014")
plot(theta14, main="5km, 3, actual 2014")

plot(psi.pred.raster15, main="5km, 3, predictions 2015")
plot(theta15, main="5km, 3, actual 2015")

plot(psi.pred.raster16, main="5km, 3, predictions 2016")
plot(theta16, main="5km, 3, actual 2016")

plot(psi.pred.raster17, main="5km, 3, predictions 2017")
plot(theta17, main="5km, 3, actual 2017")

plot(psi.pred.raster18, main="5km, 3, predictions 2018")
plot(theta18, main="5km, 3, actual 2018")

plot(psi.pred.raster19, main="5km, 3, predictions 2019")
plot(theta19, main="5km, 3, actual 2019")

#### New section on outputting results to save

toSave <- list(model = result
               ,tpr2011 = sum(preds1$theta==1 & values(theta11)==1, na.rm=TRUE)/sum(values(theta11)==1, na.rm=TRUE)
               ,fpr2011 = sum(preds1$theta==1 & values(theta11)==0, na.rm=TRUE)/sum(values(theta11)==0, na.rm=TRUE)
               ,tacc2011 = mean(preds1$theta==values(theta11), na.rm=TRUE)
               ,tpr2012 = sum(preds2$theta==1 & values(theta12)==1, na.rm=TRUE)/sum(values(theta12)==1, na.rm=TRUE)
               ,fpr2012 = sum(preds2$theta==1 & values(theta12)==0, na.rm=TRUE)/sum(values(theta12)==0, na.rm=TRUE)
               ,tacc2012 = mean(preds2$theta==values(theta12), na.rm=TRUE)
               ,tpr2013 = sum(preds3$theta==1 & values(theta13)==1, na.rm=TRUE)/sum(values(theta13)==1, na.rm=TRUE)
               ,fpr2013 = sum(preds3$theta==1 & values(theta13)==0, na.rm=TRUE)/sum(values(theta13)==0, na.rm=TRUE)
               ,tacc2013 = mean(preds3$theta==values(theta13), na.rm=TRUE)
               ,tpr2014 = sum(preds4$theta==1 & values(theta14)==1, na.rm=TRUE)/sum(values(theta14)==1, na.rm=TRUE)
               ,fpr2014 = sum(preds4$theta==1 & values(theta14)==0, na.rm=TRUE)/sum(values(theta14)==0, na.rm=TRUE)
               ,tacc2014 = mean(preds4$theta==values(theta14), na.rm=TRUE)
               ,tpr2015 = sum(preds5$theta==1 & values(theta15)==1, na.rm=TRUE)/sum(values(theta15)==1, na.rm=TRUE)
               ,fpr2015 = sum(preds5$theta==1 & values(theta15)==0, na.rm=TRUE)/sum(values(theta15)==0, na.rm=TRUE)
               ,tacc2015 = mean(preds5$theta==values(theta15), na.rm=TRUE)
               ,tpr2016 = sum(preds6$theta==1 & values(theta16)==1, na.rm=TRUE)/sum(values(theta16)==1, na.rm=TRUE)
               ,fpr2016 = sum(preds6$theta==1 & values(theta16)==0, na.rm=TRUE)/sum(values(theta16)==0, na.rm=TRUE)
               ,tacc2016 = mean(preds6$theta==values(theta16), na.rm=TRUE)
               ,tpr2017 = sum(preds7$theta==1 & values(theta17)==1, na.rm=TRUE)/sum(values(theta17)==1, na.rm=TRUE)
               ,fpr2017 = sum(preds7$theta==1 & values(theta17)==0, na.rm=TRUE)/sum(values(theta17)==0, na.rm=TRUE)
               ,tacc2017 = mean(preds7$theta==values(theta17), na.rm=TRUE)
               ,tpr2018 = sum(preds8$theta==1 & values(theta18)==1, na.rm=TRUE)/sum(values(theta18)==1, na.rm=TRUE)
               ,fpr2018 = sum(preds8$theta==1 & values(theta18)==0, na.rm=TRUE)/sum(values(theta18)==0, na.rm=TRUE)
               ,tacc2018 = mean(preds8$theta==values(theta18), na.rm=TRUE)
               ,tpr2019 = sum(preds9$theta==1 & values(theta19)==1, na.rm=TRUE)/sum(values(theta19)==1, na.rm=TRUE)
               ,fpr2019 = sum(preds9$theta==1 & values(theta19)==0, na.rm=TRUE)/sum(values(theta19)==0, na.rm=TRUE)
               ,tacc2019 = mean(preds9$theta==values(theta19), na.rm=TRUE)
)

#### NOTE: in this sample 'preds' is not year specific--make sure the right 'preds' are used when 
#### predictions are extended to multiple years

saveRDS(toSave, here("Data/Model Predictions/2001-2010_5km_5km_3_INLA.RDS"))

