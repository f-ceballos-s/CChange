library(raster)
library(rgdal)
library(rgeos)
library(parallel)
library(doSNOW)
library(snow)
library(reshape2)
library(sf)
library(tidyverse)
library(stringr)
library(velox)

## Functions

my.mean = function(x){ mean(x,na.rm=TRUE) }
my.sd = function(x){ sd(x,na.rm=TRUE) }

##

numParallelCores <- max(1, detectCores()-1)
cl <- makeCluster(rep("localhost", numParallelCores), type = "SOCK")
registerDoSNOW(cl)

setwd("C:/Users/fc3/Box Sync/Ag_Sys")

## Read Municipality Shapefile
## Read yield data
## Readg coffee producing municipalities

poly <- readOGR("C:/Users/fc3/Box Sync/Ag_Sys/Maps/MGN_MPIO_POLITICO.shp")
poly@data$MPIO_CCDGO = gsub("(^|[^0-9])0+", "",as.character(poly@data$MPIO_CCDGO))
mpio = poly@data$MPIO_CCDGO
mp_cf = read.csv("mp_cf.csv")
dem = raster("C:/Users/fc3/Box Sync/Ag_Sys/Maps/dem.tif")
yield = read.csv("yield.csv")

# Extract temperature and precipitation 
### units of temperature: K/10
### units of precipitation: kg/m^2 = mm/m^2 (monthly)

nm.dt = c("prec","tmean","tmax","tmin")

s = list()
for (i in 1:length(nm.dt)) {
  grids <- list.files(paste0("C:/Users/fc3/Box Sync/Ag_Sys/Maps/",nm.dt[i]), pattern = "*.tif$")
  s[[i]] <- stack(paste0("C:/Users/fc3/Box Sync/Ag_Sys/Maps/",nm.dt[i],"/", grids))
}

pb <- txtProgressBar(min = 0, max = length(nm.dt)*14, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
x <- foreach(i=1:4, .combine='cbind', .packages = c("raster","velox","Matrix")) %:%
  foreach(j=1:14, .combine='cbind',
          .options.snow = opts) %dopar% {
            ras_crop = crop(s[[i]]@layers[[j]],poly)
            vx.ras = velox(ras_crop)
            ex = vx.ras$extract(poly, fun=my.mean,small=T)
            return(ex)
          }
close(pb)

w <- foreach(i=1:4, .combine='c', .packages = c("raster")) %:%
  foreach(j=1:14, .combine='c',
          .options.snow = opts) %dopar% {
            name = gsub(".*CHELSA_","",s[[i]]@layers[[j]]@file@name)
            name = gsub("_V1.*","",name)
            return(name)
          }

## Create separate dataframes
## Convert temperature to C 
## wide to long 

ind_df = list()
for (i in 0:3) {
  ind_df[[i+1]] = as.data.frame(cbind(mpio,x[,((14*i)+1):((14*i)+14)]))
  colnames(ind_df[[i+1]]) = c("id",w[((14*i)+1):((14*i)+14)])
  ind_df[[i+1]][,2:15] = apply(ind_df[[i+1]][,2:15], 2, as.numeric)
}

for (i in 2:4) {
  for (j in 1:14) {
    dd = function(n){(n/10)-273}
    ind_df[[i]][,1+j] = dd(ind_df[[i]][,1+j])
  }
}

for (i in 1:4) {
  for (j in 0:1) {
    d = as.data.frame(ind_df[[i]])[,c(1,seq(2+j,14+j,2))]
    a <- melt(d,
                          id.vars="id",
                          measure.vars=colnames(d)[2:8],
                          variable.name="year",
                          value.name=paste0(nm.dt[i],"_",j)
    )
    if (i == 1 & j == 0){
      b = a
    }
    if (i != 1 | j != 0){
      b = cbind(b,a[,3,drop=FALSE])
    }
  }
}
b$year = gsub(".*prec_","",b$year)
b$year = gsub("_03.*","",b$year)

## Evapotranspiration

## Mean EVPT
bn = c(75,81,87,93,99,105,111,117)
pet_00 = brick("C:/Users/fc3/Box Sync/Ag_Sys/Maps/pet/cru_ts4.04.2001.2010.pet.dat.nc")
pet_00 = pet_00[[bn]]

bn = c(3,9,15,21,27,33)
pet_10 = brick("C:/Users/fc3/Box Sync/Ag_Sys/Maps/pet/cru_ts4.04.2011.2019.pet.dat.nc")
pet_10 = pet_10[[bn]]

pet = stack(pet_00,pet_10)

pb <- txtProgressBar(min = 0, max = 14, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

x = as.data.frame(foreach(j=1:length(pet@layers), .combine='cbind',
                           .options.snow = opts, .packages = c("raster","velox","Matrix")) %dopar% {
                             ras_crop = crop(pet@layers[[j]],poly)
                             vx.ras = velox(ras_crop)
                             ex = vx.ras$extract(poly, fun=my.mean,small=T)
                             return(ex)
                           })

w = foreach(j=1:14, .combine='c',
        .options.snow = opts, .packages = "raster") %dopar% {
          name = names(pet@layers[[j]])
          return(name)
        }

colnames(x) = w
x$id = poly@data$MPIO_CCDGO

z = reshape(x, direction='long', 
            varying=list(colnames(x)[seq(1,14,2)],colnames(x)[seq(2,14,2)]), 
            timevar="year",
            times=2007:2013,
            v.names=c("pet_mar","pet_aug"),
            idvar='id')

## Altitude

vx.ras = velox(dem)
alt = vx.ras$extract(poly, fun=my.mean,small=T)
alt_sd = vx.ras$extract(poly, fun=my.sd,small=T)
alt[which(is.na(alt))] = 0 # San Andres y Providencia
altitude = as.data.frame(cbind(poly@data$MPIO_CCDGO,alt,alt_sd))
colnames(altitude) = c("id","alt_mn","alt_sd")

## Merge all

df = Reduce(function(x, y) merge(x, y, by=c("id","year"), all=TRUE), list(yield,b,z))
df = left_join(df,altitude,by="id",all.y=F)
df = df[which(df$id %in% mp_cf$mp_cf),]
colnames(df)[6:ncol(df)] = c("prec_mar","prec_aug","tmean_mar","tmean_aug","tmax_mar","tmax_aug","tmin_mar","tmin_aug","pet_mar","pet_aug","alt_mn","alt_sd")

## Create altitude dummy
## Create split data sets

df = apply(df, 2, as.numeric)
df = as.data.frame(cbind(df,ifelse(df[,"alt_mn"] > mean(df[,"alt_mn"],na.rm = T),1,0)))
colnames(df)[ncol(df)] = "alt_dm"
hi_df = df[which(df[,"alt_dm"]==1),]
lo_df = df[which(df[,"alt_dm"]==0),]
df$dpto = gsub('.{3}$', '', df$id)

write.csv(df,"df.csv",row.names = F)

stopCluster(cl)
