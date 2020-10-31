library(raster)
library(rgdal)
library(rgeos)
library(parallel)
library(doSNOW)
library(snow)
library(velox)
library(tidyverse)

numParallelCores <- max(1, detectCores()-1)
cl <- makeCluster(rep("localhost", numParallelCores), type = "SOCK")
registerDoSNOW(cl)

my.mean = function(x){ mean(x,na.rm=TRUE) }
my.sd = function(x){ sd(x,na.rm=TRUE) }
t_10 = function(x){x/10}

setwd("C:/Users/fc3/Box Sync/Ag_Sys")

# Functions

## Read Municipality Shapefile
## Read yield data
## Readg coffee producing municipalities

poly <- readOGR("C:/Users/fc3/Box Sync/Ag_Sys/Maps/MGN_MPIO_POLITICO.shp")
poly@data$MPIO_CCDGO = gsub("(^|[^0-9])0+", "",as.character(poly@data$MPIO_CCDGO))
dem = raster("C:/Users/fc3/Box Sync/Ag_Sys/Maps/dem.tif")

# Extract tmean and prec
### units of temperature: K/10
### units of precipitation: kg/m^2 = mm/m^2 (monthly)

setwd("C:/Users/fc3/Box Sync/Ag_Sys/Maps/proj2")

nm.dt = c("BCC-CSM1","CNRM-CM5","CanESM2",
          "IPSL-CM5A-LR","MIROC-ESM","MIROC5","MRI-CGCM3")

s = list()
pb <- txtProgressBar(min = 0, max = length(nm.dt), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
x <- as.data.frame(foreach(i=1:length(nm.dt), .combine='cbind', .packages = c("raster","velox", "Matrix")) %:%
                     foreach(j=1:4, .combine='cbind',
                             .options.snow = opts) %dopar% {
                               grids <- list.files(paste0("C:/Users/fc3/Box Sync/Ag_Sys/Maps/proj2/",nm.dt[i]), pattern = "*.tif$")
                               s[[i]] <- stack(paste0("C:/Users/fc3/Box Sync/Ag_Sys/Maps/proj2/",nm.dt[i],"/", grids))
                               ras_crop = crop(s[[i]]@layers[[j]],poly)
                               ras_crop[ras_crop < 0 | ras_crop > 1000] = NA
                               # ras_crop[ras_crop > 500] = 500
                               vx.ras = velox(ras_crop)
                               ex = vx.ras$extract(poly, fun=my.mean,small=T)
                               return(ex)
                             })

close(pb)
x = cbind(id = poly@data$MPIO_CCDGO,x)
w = reshape(x, direction='long', 
            varying=list(colnames(x)[seq(2,29,4)],colnames(x)[seq(3,29,4)],colnames(x)[seq(4,29,4)],colnames(x)[seq(5,29,4)]), 
            timevar="model",
            times=c("BCC-CSM1","CNRM-CM5","CanESM2","IPSL-CM5A-LR","MIROC-ESM","MIROC5","MRI-CGCM3"),
            v.names=c("prec_mar", "prec_aug","tmean_mar","tmean_aug"),
            idvar='id')
colnames(w)[1] = "id"
w$id = as.character(w$id)

w[,c("tmean_mar","tmean_aug")] = sapply(w[,c("tmean_mar","tmean_aug")], t_10)

# Altitude

vx.ras = velox(dem)
alt = vx.ras$extract(poly, fun=my.mean,small=T)
alt[which(is.na(alt))] = 0 # San Andres y Providencia
altitude = as.data.frame(cbind(id = poly@data$MPIO_CCDGO,alt))
colnames(altitude) = c("id","alt_mn")


df = left_join(w,altitude,by="id",all=F)

setwd("C:/Users/fc3/Box Sync/Ag_sys")

write.csv(df,"proj2_df.csv", row.names = F)

stopCluster(cl)
