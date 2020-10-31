setwd("C:/Users/fc3/Box/CChange (fc3@illinois.edu)")

library(tidyverse)
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(psych)
library(grid)
library(Gmisc)
library(foreach)
library(rgdal)
library(rgeos)
library(viridis)
library(sf)
library(gridExtra)
library(png)

bar_gph = list()
ggproj = list()
chg.ls = list()

rcp = c("proj","proj2","proj3")
l = 3
chg = matrix(0,3,3)

setwd("C:/Users/fc3/Box/Ag_sys")

df = read.csv("df.csv")
pr_df = read.csv(paste0(rcp[l],"_df.csv"))
poly <- readOGR("C:/Users/fc3/Box/Ag_Sys/Maps/MGN_MPIO_POLITICO.shp")
poly@data$MPIO_CCDGO = gsub("(^|[^0-9])0+", "",as.character(poly@data$MPIO_CCDGO))

################################
########### Functions###########
################################

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


y.s = function(tmean_mar,tmean_aug,prec_mar,prec_aug,alt_mn){
  y = .8902765*tmean_mar - .0178434*tmean_mar^2 + .2138056*tmean_aug - .0179668*tmean_aug^2 + 
    .010204*prec_mar + 0.00000896*prec_mar^2 - .0056075*prec_aug - 0.00000133*prec_aug^2 - 
    .0005797*tmean_mar*prec_mar + 0.000327*tmean_aug*prec_aug - .001743*alt_mn - 4.144
  return(y)}

confidence_interval <- function(vector,weight, interval) {
  vec_sd <- sd(vector, na.rm = T)
  n <- length(!is.na(vector))
  vec_mean <- weighted.mean(vector,weight, na.rm = T)
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  result <- c("mean" = vec_mean,"lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

convert_func<-function(x){  as.numeric(as.character(x))}

cut_with_nas   <- function( x, breaks, labels=NULL, .missing="Unknown" ) {
  y <- cut(x, breaks, labels) #, include.lowest = T, right=F)
  y <- addNA(y)
  levels(y)[is.na(levels(y))] <- .missing
  return( y )
}

pd <- position_dodge(width = 0.4)

################################
######### Projections###########
################################

df$as_cafe[df$as_cafe == 0] = 0.1

for (i in 1:nrow(df)) {
  df$v_cafe[i] = df$p_cafe[i]/df$as_cafe[i]
}

df = df %>%
  dplyr::group_by(id)%>%
  dplyr::mutate(lag.value = dplyr::lag(v_cafe, n = 1, default = NA))

## Sys GMM
df$y_hat.s = y.s(df$tmean_mar,df$tmean_aug,df$prec_mar,df$prec_aug,df$alt_mn)
df$res.s = df$v_cafe - df$y_hat.s
c_i = aggregate(df$res.s, list(df$id), mean)


dif_min = list()
sys_max = list()
sys_min = list()
sys_max = list()
files <- list.files("C:/Users/fc3/Box/CChange (fc3@illinois.edu)/Maps/proj/proj_df", pattern = "*.csv$")
temp = c("tmax","tmin")
for (i in 1:length(files)) {
  a <- read.csv(paste0("C:/Users/fc3/Box/CChange (fc3@illinois.edu)/Maps/proj/proj_df/", files[i]))
  b = df[which(df$year==2013),c("id","alt_mn","alt_dm")]
  a = left_join(a, b, by = c("MPIO_CCDGO" = "id"))
  a = left_join(a, c_i, by = c("MPIO_CCDGO" = "Group.1"))
  a$tmean_mar = (a$tmax_mar+a$tmin_mar)/2
  a$tmean_aug = (a$tmax_aug+a$tmin_aug)/2
  attach(a)
  a$proj = y.s(tmean_mar,tmean_aug,prec_mar,prec_aug,alt_mn)
  a$proj = a$proj
  
  # for (k in 1:length(alt_dm)) {
  #   if(tmean_mar[k] > 28.515320 | tmean_mar[k] < 10.11 |
  #      tmean_aug[k] > 27.917779 | tmean_aug[k] < 9.98 |
  #      prec_mar[k] > 481.951351 | prec_mar[k] < 7.728571 |
  #      prec_aug[k] > 491.662162 | prec_aug[k] < 1.681818){a[k,c("proj")] <- NA}
  # }
  sys_max[[i]] = a
}

###############################################
########## Descriptive statistics##############
###############################################

## Present

int.vars = c("v_cafe","as_cafe","p_cafe","prec_mar","tmean_mar","prec_aug","tmean_aug","alt_mn","alt_dm")
dsc_df = df[,int.vars]
descr = describeBy(dsc_df[,1:7],factor(dsc_df$alt_dm),mat = T)
dsc_stat_pres = cbind(descr[which(descr$group1 == 0),c("mean","sd","min","max")],descr[which(descr$group1 == 1),c("mean","sd","min","max")])
colnames(dsc_stat_pres)[((ncol(dsc_stat_pres)/2+1)):ncol(dsc_stat_pres)] = paste0(colnames(dsc_stat_pres)[((ncol(dsc_stat_pres)/2+1)):ncol(dsc_stat_pres)],"_2")


int.vars = c("prec_mar","prec_aug","tmax_mar","tmax_aug","tmin_mar","tmin_aug","alt_mn","alt_dm")
alt_dm = sys_max[[1]]$alt_dm
dsc_stat_fut = matrix(0,length(int.vars),8)
for (i in 1:length(int.vars)) {
  ## Max
  a = lapply(sys_max, "[[", int.vars[i])
  maxd = do.call("cbind", a)
  dsc_stat_fut[i,1:4] = c(mean(rowSums(maxd[which(alt_dm == 0),])/8),sd(rowSums(maxd[which(alt_dm == 0),])/8),min(rowSums(maxd[which(alt_dm == 0),])/8),max(rowSums(maxd[which(alt_dm == 0),]))/8)
  dsc_stat_fut[i,5:8] = c(mean(rowSums(maxd[which(alt_dm == 1),])/8),sd(rowSums(maxd[which(alt_dm == 1),])/8),min(rowSums(maxd[which(alt_dm == 1),])/8),max(rowSums(maxd[which(alt_dm == 1),]))/8)
}
dsc_stat_fut = as.data.frame(dsc_stat_fut)
rownames(dsc_stat_fut) = int.vars
colnames(dsc_stat_fut) = colnames(dsc_stat_pres)
dsc_stat = rbind(dsc_stat_pres,dsc_stat_fut)

dsc_stat = apply(dsc_stat, 2, round,2)

################################
########### line plots##########
################################

cl.nm = c("prec_mar","prec_aug","tmean_mar","tmean_aug")

ln.gph = list()
m_df = reshape(as.data.frame(df[,c("id","year","prec_mar","prec_aug","tmean_mar","tmean_aug")]), direction='long', 
               varying= list(c("prec_mar","prec_aug"),c("tmean_mar","tmean_aug")), 
               timevar="month",
               times=c("March","August"),
               v.names= c("Precipitation","Temperature"),
               idvar=c("id","year"))
a = data_summary(data=m_df,varname="Temperature",groupnames=c("year","month"))
c = data_summary(data=m_df,varname="Precipitation",groupnames=c("year","month"))
b = data_summary(data=df,varname="p_cafe",groupnames=c("year"))

nh = c("a","c")

## Create graphs
ln.gph[[1]] = ggplot(a, aes(x=factor(year), y=Temperature, group=month, color=month)) +
  geom_errorbar(aes(ymin=Temperature-sd, ymax=Temperature+sd), width=.1, position = pd) +
  geom_line(aes(x=factor(year), y=Temperature), position = pd) +
  geom_point(aes(x = factor(year), y = Temperature, shape=factor(month)), position = pd)+
  scale_color_manual(name="Series",
                     labels = c("August temperature", "March temperature"), 
                     values = c("gray60", "gray8")) +
  scale_shape_manual(name="Series",
                     labels = c("August temperature", "March temperature"), 
                     values = c(13,15)) +
  ylab("Temperature (°C)")+
  xlab("Year")+
  guides(linetype=FALSE,color=FALSE)+
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.text=element_text(size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 16),
        panel.background = element_rect(fill = "#FFFFFF", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray70"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray70"))

ln.gph[[2]] = ggplot(c, aes(x=factor(year), y=Precipitation, group=month, color=month)) +
  geom_errorbar(aes(ymin=Precipitation-sd, ymax=Precipitation+sd), width=.1, position = pd) +
  geom_line(aes(x=factor(year), y=Precipitation), position = pd) +
  geom_point(aes(x = factor(year), y = Precipitation, shape=factor(month)), position = pd)+
  scale_color_manual(name="Series",
                     labels = c("August precipitation", "March precipitation"), 
                     values = c("gray60", "gray8")) +
  scale_shape_manual(name="Series",
                     labels = c("August precipitation", "March precipitation"), 
                     values = c(13,15)) +
  ylab("Precipiation (mm/month)")+
  xlab("Year")+
  guides(linetype=FALSE,color=FALSE)+
  theme(legend.position="bottom",
        legend.direction="horizontal",
        legend.text=element_text(size=16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 16),
        panel.background = element_rect(fill = "#FFFFFF", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray70"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray70"))

ln.gph[[3]] = ggplot(b, aes(x=year, y=p_cafe)) +
  geom_line(aes(x=year, y=p_cafe), position = pd) +
  geom_point(aes(x = year, y = p_cafe), position = pd)+
  ylab("National coffee yield (tons)")+
  xlab("Year")+
  guides(linetype=FALSE,color=FALSE)+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        panel.background = element_rect(fill = "#FFFFFF", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray70"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray70"))

ggline = as_ggplot(grid.arrange(arrangeGrob(ln.gph[[1]], ln.gph[[2]], ncol = 2),                             # First row with one plot spaning over 2 columns
                                ln.gph[[3]], # Second row with 2 plots in 2 different columns
                                nrow = 2))
ggsave("ggline2.png", width = 40, height = 20, units = "cm")

## Temperature in March

x <- 10:30
dat <- data.frame(x,
                  ys1= y.s(x,median(df$tmean_aug[which(df$alt_dm ==1)]),median(df$prec_mar[which(df$alt_dm ==1)]),median(df$prec_aug[which(df$alt_dm ==1)]),median(df$alt_mn[which(df$alt_dm ==1)],na.rm=T)),
                  ys2= y.s(x,median(df$tmean_aug[which(df$alt_dm ==0)]),median(df$prec_mar[which(df$alt_dm ==0)]),median(df$prec_aug[which(df$alt_dm ==0)]),median(df$alt_mn[which(df$alt_dm ==0)],na.rm=T))
)
dat$upr1 = dat$ys1+0.52
dat$lwr1 = dat$ys1-0.52
dat$upr2 = dat$ys2+0.36
dat$lwr2 = dat$ys2-0.36
ggtemp_m = ggplot(dat, aes(x)) + 
  geom_point(aes(y=ys1, colour = "System GMM (high altitude)"),shape=15, size = 2)+
  geom_line(aes(y=ys1))+
  geom_ribbon(data=dat,aes(ymin=lwr1,ymax=upr1),alpha=0.15)+
  geom_point(aes(y=ys2, colour = "System GMM (low altitude)"),shape=17, size = 2)+
  geom_line(aes(y=ys2))+
  geom_ribbon(data=dat,aes(ymin=lwr2,ymax=upr2),alpha=0.25)+
  geom_point(aes(x=median(df$tmean_mar[which(df$alt_dm ==1)]), y=y.s(median(df$tmean_mar[which(df$alt_dm ==1)]),median(df$tmean_aug[which(df$alt_dm ==1)]),median(df$prec_mar[which(df$alt_dm ==1)]),median(df$prec_aug[which(df$alt_dm ==1)]),median(df$alt_mn[which(df$alt_dm ==1)],na.rm=T))), colour="black", shape=8, size = 8) +
  annotate(geom = "text", x = 14.5, y = 1.83, label = "(17.17,1.88)", color = "black",angle = 0, size = 3) +
  geom_point(aes(x=median(df$tmean_mar[which(df$alt_dm ==0)]), y=y.s(median(df$tmean_mar[which(df$alt_dm ==0)]),median(df$tmean_aug[which(df$alt_dm ==0)]),median(df$prec_mar[which(df$alt_dm ==0)]),median(df$prec_aug[which(df$alt_dm ==0)]),median(df$alt_mn[which(df$alt_dm ==0)],na.rm=T))), colour="black", shape=8, size = 8) +
  annotate(geom = "text", x = 22, y = 1.85, label = "(22.05,1.68)", color = "black",angle = 0, size = 3) +
  geom_vline(xintercept = min(df$tmean_mar),linetype="dotted")+
  geom_vline(xintercept = max(df$tmean_mar),linetype="dotted")+ 
  labs(x ="Degree Celsius (°C)", y = "Coffee productivity (1000 kg./ha.)")+
  scale_color_manual(name = "Altitude set",
                     values = c("black","black"))+
  guides(colour = guide_legend(override.aes=list(shape = c(15,17))))+
  ylim(c(-1.5,3))+  
  theme(legend.text=element_text(size=20))+
  theme_bw()


## Precipitation in March

x <- seq(0,550,20)
dat <- data.frame(x,
                  ys1= y.s(median(df$tmean_mar[which(df$alt_dm ==1)]),median(df$tmean_aug[which(df$alt_dm ==1)]),x,median(df$prec_aug[which(df$alt_dm ==1)]),median(df$alt_mn[which(df$alt_dm ==1)],na.rm=T)),
                  ys2= y.s(median(df$tmean_mar[which(df$alt_dm ==0)]),median(df$tmean_aug[which(df$alt_dm ==0)]),x,median(df$prec_aug[which(df$alt_dm ==0)]),median(df$alt_mn[which(df$alt_dm ==0)],na.rm=T))
)
dat$upr1 = dat$ys1+0.52
dat$lwr1 = dat$ys1-0.52
dat$upr2 = dat$ys2+0.36
dat$lwr2 = dat$ys2-0.36
ggprec_m = ggplot(dat, aes(x)) + 
  geom_point(aes(y=ys1, colour = "System GMM (high altitude)"),shape=15, size = 2)+
  geom_line(aes(y=ys1))+
  geom_ribbon(data=dat,aes(ymin=lwr1,ymax=upr1),alpha=0.15)+
  geom_point(aes(y=ys2, colour = "System GMM (low altitude)"),shape=17, size = 2)+
  geom_line(aes(y=ys2))+
  geom_ribbon(data=dat,aes(ymin=lwr2,ymax=upr2),alpha=0.15)+
  geom_vline(xintercept = min(df$prec_mar),linetype="dotted")+
  geom_vline(xintercept = max(df$prec_mar),linetype="dotted")+
  geom_point(aes(x=median(df$prec_mar[which(df$alt_dm ==1)]), y=y.s(median(df$tmean_mar[which(df$alt_dm ==1)]),median(df$tmean_aug[which(df$alt_dm ==1)]),median(df$prec_mar[which(df$alt_dm ==1)]),median(df$prec_aug[which(df$alt_dm ==1)]),median(df$alt_mn[which(df$alt_dm ==1)],na.rm=T))), colour="black", shape=8, size = 6) +
  annotate(geom = "text", x = 110, y = 2, label = "(119.21,1.88)", color = "black",angle = 0, size = 3) +
  geom_point(aes(x=median(df$prec_mar[which(df$alt_dm ==0)]), y=y.s(median(df$tmean_mar[which(df$alt_dm ==0)]),median(df$tmean_aug[which(df$alt_dm ==0)]),median(df$prec_mar[which(df$alt_dm ==0)]),median(df$prec_aug[which(df$alt_dm ==0)]),median(df$alt_mn[which(df$alt_dm ==0)],na.rm=T))), colour="black", shape=8, size = 6) +
  annotate(geom = "text", x = 100, y = 1.6, label = "(108.69,1.68)", color = "black",angle = 0, size = 3) +
  labs(x ="mm.", y = "Coffee productivity (1000 kg./ha.)")+
  scale_color_manual(name = "Altitude set",
                     values = c("black","black"))+
  guides(colour = guide_legend(override.aes=list(shape = c(15,17))))+
  ylim(c(-1,3))+  
  theme(legend.text=element_text(size=20))+
  theme_bw()

ggmarg = ggarrange(ggtemp_m, ggprec_m,
                   ncol = 2, nrow = 1, 
                   legend="bottom", common.legend = TRUE,
                   labels = c("A","B"))

ggsave("ggmarg2.png", width = 40, height = 20, units = "cm")

#################################################################
########## Projected climate: graphs and map#####################
#################################################################

## projections

pr_df = left_join(pr_df,df[which(df$year == 2013),c("id","as_cafe","v_cafe")],by="id")
pr_df = pr_df[complete.cases(pr_df[,"as_cafe"]),]
pr_df$alt_dm = ifelse(pr_df$alt_mn > median(pr_df$alt_mn),1,0)
pr_df$proj = y.s(pr_df$tmean_mar,pr_df$tmean_aug,pr_df$prec_mar,pr_df$prec_aug,pr_df$alt_mn)
pr_df$proj[pr_df$proj < 0] = 0

nm.dt = c("BCC-CSM1","CNRM-CM5","CanESM2",
          "IPSL-CM5A-LR","MIROC-ESM","MIROC5","MRI-CGCM3")
x = foreach(i = 1:length(nm.dt), .combine = "rbind")%do%{
  confidence_interval(pr_df$proj[which(pr_df$model == nm.dt[i])],pr_df$as_cafe[which(pr_df$model == nm.dt[i])],.95)
}

y = foreach(j=0:1, .combine = "rbind")%:%
  foreach(i=1:length(nm.dt), .combine = "rbind")%do%{
    confidence_interval(pr_df$proj[which(pr_df$model == nm.dt[i] & pr_df$alt_dm == j)],pr_df$as_cafe[which(pr_df$model == nm.dt[i] & pr_df$alt_dm == j)],.95)
  }

z = as.data.frame(rbind(x,y))
z$model = rep(c("BCC-CSM1","CNRM-CM5","CanESM2","IPSL-CM5A-LR","MIROC-ESM","MIROC5","MRI-CGCM3"),3)
z$alt = c(rep("all",7),rep("low",7),rep("high",7))


pd <- position_dodge(width = 0.2)
bar_gph[[l]] = ggplot(z, aes(x=model, y=mean, group=alt, color=alt, shape = alt)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.15, position = pd) +
  geom_point(aes(x = model, y = mean), position = pd)+
  geom_hline(yintercept = mean(df$v_cafe,na.rm=T, linetype = "Mean productivity"))+
  scale_color_manual(name = "Altitude \ngroup",
                     labels = c("Whole sample","High altitude","Low altitude"),
                     values = c("all"="brown1","low"="darkgoldenrod1","high"="dodgerblue1"))+
  scale_shape_manual(name = "Altitude \ngroup",
                     labels = c("Whole sample","High altitude","Low altitude"),
                     values = c("all"=13,"low"=15,"high"=17)) +
  xlab("Global Circulation Models (GCM)") +
  ylab("Coffee productivity (1000 kg./ha)")+
  # ylim(c(0,3)) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        panel.background = element_rect(fill = "#FFFFFF", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray70"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray70"))

chg[1,] = as.vector(confidence_interval(pr_df[,"proj"],pr_df[,"as_cafe"],.95))
chg[2,] = as.vector(confidence_interval(pr_df[which( pr_df$alt_dm == 0),"proj"],pr_df[which( pr_df$alt_dm == 0),"as_cafe"],.95))
chg[3,] = as.vector(confidence_interval(pr_df[which( pr_df$alt_dm == 1),"proj"],pr_df[which( pr_df$alt_dm == 1),"as_cafe"],.95))

chg.ls[[l]] = as.data.frame(chg)

if(l == 3){chg_df = bind_rows(chg.ls[[3]],chg.ls[[2]],chg.ls[[1]])

colnames(chg_df) = c("mean","lower","upper")
chg_df$alt = rep(c("all","low","high"),3)
chg_df$rcp = c(rep("rcp6.0",3),rep("rcp2.6",3),rep("rcp4.5",3))

}

bar_hor = ggplot(chg_df, aes(x=factor(rcp), y=mean, group=factor(alt), color=factor(alt), shape = factor(alt))) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.15, position = pd) +
  geom_point(aes(x = factor(rcp), y = mean), position = pd)+
  geom_hline(yintercept = mean(df$v_cafe,na.rm=T, linetype = "Mean productivity"), color = "black")+
  scale_color_manual(name = "Altitude \ngroup",
                     labels = c("Whole sample","Low altitude","High altitude"),
                     values = c("all"="brown1","low"="darkgoldenrod1","high"="dodgerblue1"))+
  scale_shape_manual(name = "Altitude \ngroup",
                     labels = c("Whole sample","High altitude","Low altitude"),
                     values = c("all"=13,"low"=15,"high"=17)) +
  xlab("Representative Concentration Pathway (RCP)") +
  ylab("Coffee productivity (1000 kg./ha)")+
  # ylim(c(0,3)) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        panel.background = element_rect(fill = "#FFFFFF", colour = "black",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "gray70"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray70"))

if(l == 3){ggbar = ggarrange(bar_gph[[2]], bar_gph[[1]], bar_gph[[3]], bar_hor,
                   ncol = 2, nrow = 2,
                   legend="right", common.legend = TRUE,
                   labels = c("A","B","C","D"))

ggsave("ggbar.png", width = 50, height = 40, units = "cm")}

############################################################# Map############################################################################
# 
# pr_df$chg = pr_df$proj - pr_df$v_cafe
# pr_df$id = as.character(pr_df$id)
# 
# v_df = ggplot2::fortify(poly, region = "DPTOMPIO", group = "DPTOMPIO")
# 
# v.mn = v_df
# v_df = left_join(v_df, pr_df, by = c("id"))
# v_df$chg_v = ifelse(v_df$chg < 0,"darkblue","yellow")
# ggproj[[l]] = ggplot()+
#   geom_polygon(data = v_df,
#                aes(long, lat, group=group, fill = chg_v),
#                color="black", size=0.01)+
#   scale_fill_identity(name = "",
#                       labels = c("Expected \ndecrease","Expected \nincrease"),
#                       guide = "legend")+
#   theme_map()+
#   theme(legend.position="left",legend.text=element_text(size=12))+
#   labs(color='% change in night lights, 2012-2016')
# 
# if(l == 3){ggmap = ggarrange(ggproj[[2]], ggproj[[1]], ggproj[[3]],
#                    ncol = 3, nrow = 1,
#                    legend="bottom", common.legend = TRUE,
#                    labels = c("A","B","C"))
# 
# ggsave("ggmap.png", width = 40, height = 20, units = "cm")}

